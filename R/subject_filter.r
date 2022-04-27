# Functions for generating list of subjects to keep based on their mutation type signatures
library(NMF)
library(smaug)
library(tidyverse)

get_mct_b <- function(bins){
  out <- bins %>%
    group_by(CHR, BIN, Motif) %>%
    summarise(nMotifs=sum(nMotifs))
  return(out)
}

sub_mer <- function(motif, n){
  lMotif <- (nchar(motif) - 2) / 2
  center <- ceiling(lMotif / 2)
  SEQA <- substr(motif, center - n, center + n)
  SEQB <- substr(motif, (center * 3) - n, (center * 3) + n)
  return(paste0(SEQA, "(", SEQB, ")"))
}

get_ind_data <- function(pedfile, phenofile, qplotfile, contamfile, vcfastfile, pcfile){
  options(readr.show_progress = FALSE)

  # get ids in ped file
  ped <- read_tsv(pedfile)
  names(ped)[1] <- "ID"

  pheno <- read_tsv(phenofile)
  pheno <- pheno %>%
    dplyr::select(ID=SEQID, Sex, Case.Control, Study, DROP) %>%
    filter(!(is.na(ID)))
  # RESULT : ID, Sex, Case.Control, Study, DROP

  qplotdat <- read_tsv(qplotfile)
  names(qplotdat)[1] <- "ID"

  contam <- read_table(contamfile)
  names(contam) <- c("ID", "FREEMIX", "CHIPMIX", "PLATE")

  vcfast <- read.table(vcfastfile,
                       header=T, stringsAsFactors=F, sep="\t", comment.char=" ")
  names(vcfast) <- c("ID", "SNPs", "Singletons", "Doubletons", "lt0.5", "gt0.5", "Ref", "Het", "Alt", "Heterozygosity")

  pcs <- read.table(pcfile, header=F, stringsAsFactors=F, skip=1)
  names(pcs) <- c("ID", paste0("PC",1:10), "Study")
  pcs <- pcs %>%
    dplyr::select(-Study)

  subjectData <- merge(pheno, pcs, by="ID") %>%
    merge(vcfast, by="ID") %>%
    merge(contam, by="ID", all.x=TRUE) %>% distinct(ID, .keep_all = TRUE)

  ped_b_ids <- subjectData[!(subjectData$ID %in% qplotdat$ID),] %>%
    dplyr::select(ID_B=ID) %>%
    mutate(ID=gsub("B", "", ID_B))

  # Subset ped2 to B ids
  qplotB <- qplotdat %>%
    filter(ID %in% ped_b_ids$ID) %>%
    mutate(ID=paste0(ID, "B"))

  qplot2 <- rbind(qplotdat, qplotB)
  subjectData <- merge(subjectData, qplot2, by="ID")
}

get_loads <- function(widedat, nmfdat){
  sigloads <- data.frame(subtype=names(widedat),
                         sig1=coef(nmfdat)[1,],
                         sig2=coef(nmfdat)[2,],
                         sig3=coef(nmfdat)[3,]) %>%
    mutate(sig1=sig1/sum(sig1), sig2=sig2/sum(sig2), sig3=sig3/sum(sig3)) %>%
    gather(sig, value, sig1:sig3)

  names(sigloads) <- c("subtype", "sig", "value")
  sigloads <- sigloads %>%
    mutate(Category=substr(subtype, 1, 5),
           Sequence=substr(subtype, 7, 14))

  return(sigloads)
}

plot_loads <- function(sigloads){
  p <- ggplot(sigloads, aes(x=Sequence, y=value, fill=sig))+
    geom_bar(stat="identity")+
    facet_grid(sig~Category, scales="free_x")+
    # geom_label_repel(data=sigloads[sigloads$sig3>0.005,], aes(x=sig2, y=sig3, label=subtype))+
    theme_bw()+
    theme(axis.text.x=element_text(angle=90, hjust=1),
          strip.text=element_text(size=16),
          legend.position="none")
  return(p)
}

plot_ind_sigs <- function(sigdat){
  p <- ggplot(sigdat,
              aes(x=ID, y=prob, group=factor(Signature), colour=factor(Signature)))+
    geom_line()+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), limits=c(0,1))+
    facet_wrap(~Study, scales="free_x", nrow=1)+
    ylab("signature contribution")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          legend.position="bottom")

  return(p)
}

plot_ind_sigs2 <- function(sigdat){
  p <- ggplot(sigdat,
              aes(x=ID, y=prob, group=factor(Signature), colour=factor(Signature)))+
    geom_line()+
    scale_x_discrete(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0), limits=c(0,1))+
    facet_wrap(~top_r, scales="free_x", nrow=1)+
    ylab("signature contribution")+
    theme_bw()+
    theme(axis.text.x=element_blank(),
          legend.position="bottom")

  return(p)
}

get_nmf_input <- function(sites, binDir, binPattern, subjectData){

  bins1Mb <- get_bins(binDir, binPattern)
  bins1Mb$CHR <- gsub("chr", "", bins1Mb$CHR)

  # Count singletons per 3-mer subtype per individual
  ind_counts <- sites %>%
    mutate(Type=gsub("cpg_", "", Category2),
           Motif=sub_mer(Motif,1)) %>%
    group_by(ID, Type, Motif) %>%
    summarise(n=n())
  # RESULT: ID, Type, Motif, n

  # count 3-mer motifs genome-wide
  mct3 <- get_mct_b(bins1Mb) %>%
    group_by(Motif) %>%
    summarise(nMotifs=sum(nMotifs))
  # RESULT : Motif, nMotif

  # get per-person rates and filter IDs
  ind_counts <- merge(ind_counts, mct3, by="Motif") %>%
    filter(ID %in% subjectData$ID) %>%
    mutate(ERV_rel_rate=n/nMotifs,
           subtype=paste0(Type, "_", Motif))

  ind_wide <- ind_counts	%>%
    dplyr::select(ID, subtype, ERV_rel_rate) %>%
    spread(subtype, ERV_rel_rate)
  ind_wide[is.na(ind_wide)] <- 0
  # RESULT: subject x rel_rates_per_3mer

  return(ind_wide)
}

nmf_analysis <- function(ind_wide, subjectData, analysisdir, outPlots = TRUE){
  ind_nmf <- nmf(as.matrix(ind_wide[,-c(1)]), 3, nrun=10)
  ind_pred <- predict(ind_nmf, "rows", prob=T)

  if(outPlots){
    sigloads <- get_loads(ind_wide[,-c(1)], ind_nmf)
    plot_loads(sigloads)
    ggsave(paste0(analysisdir, "/images/sigloads.png"), width=12, height=6)
  }

  # Add NMF basis and predicted class to data frame
  nmfdat1 <- data.frame(ID=ind_wide$ID, basis(ind_nmf), sig=ind_pred[[1]])

  nmfdat1 <- merge(nmfdat1, subjectData, by="ID") %>%
    mutate(sum=X1+X2+X3,
           sig1=X1/sum,
           sig2=X2/sum,
           sig3=X3/sum) %>%
    mutate(top_r=apply(.[,51:53], 1, function(x) names(x)[which.max(x)]))

  nmfdat1$ID <- factor(nmfdat1$ID, levels = unique(nmfdat1$ID))

  svmdat <- nmfdat1 %>%
    filter(!(is.na(PLATE))) %>%
    filter(!(is.na(CHIPMIX)))

  ind_nmf_long <- nmfdat1 %>%
    gather(Signature, prob, sig1:sig3) %>%
    mutate(Signature=as.numeric(gsub("sig", "", Signature))) %>%
    arrange(Study, sum)

  if(outPlots){
    plot_ind_sigs(ind_nmf_long)
    ggsave(paste0(analysisdir, "/images/by_ind_all.png"), width=8, height=6)

    plot_ind_sigs2(ind_nmf_long)
    ggsave(paste0(analysisdir, "/images/by_ind_all_by_sig.png"), width=8, height=6)
  }

  keep_cis <- ind_nmf_long %>%
    dplyr::select(Signature, prob) %>%
    group_by(Signature) %>%
    summarise_each(funs(mean,sd)) %>%
    mutate(conf.low=mean-1.96*sd, conf.high=mean+1.96*sd)

  nmfdat1 <- nmfdat1 %>%
    mutate(flag=ifelse(sig1 > keep_cis[1,]$conf.high, "sig1",
                       ifelse(sig3 > keep_cis[3,]$conf.high, "sig3", "sig2")))

  keep_ids <- nmfdat1 %>%
    filter(sig1 > keep_cis[1,]$conf.low & sig1 < keep_cis[1,]$conf.high) %>%
    filter(sig2 > keep_cis[2,]$conf.low & sig2 < keep_cis[2,]$conf.high) %>%
    filter(sig3 > keep_cis[3,]$conf.low & sig3 < keep_cis[3,]$conf.high) %>%
    dplyr::select(ID)

  drop_ids <- nmfdat1 %>%
    filter(!(ID %in% keep_ids$ID)) %>%
    dplyr::select(ID, flag)
  return(keep_ids)
}
