#!/usr/bin/RScript

cat("Loading packages...")
install.packages("yaml", quiet=TRUE)
require(yaml)

args <- yaml.load_file("./_config.yaml")
attach(args)

# add user libpath, if it exists (useful if running on a different node)
.libPaths(c(.libPaths(), args$libpath))

require(devtools)
install_github('theandyb/smaug', quiet=TRUE)
install_github('slowkow/ggrepel', quiet=TRUE)

gh_packages <- c("smaug", "ggrepel")
invisible(sapply(gh_packages, function(x)
	suppressMessages(require(x, character.only = TRUE))))

# install/load CRAN packages
packages <- c("tidyverse", "broom", "RColorBrewer", "MASS", "boot", "speedglm",
	"psych", "lmtest", "fmsb", "hexbin", "cowplot", "grid", "gtable", "gridExtra",
	"yaml", "openxlsx", "Biostrings", "svglite", "NMF", "emdbook")
invisible(sapply(packages, function(x)
	suppressMessages(usePackage(x))))

# Define additional variables for cleaner strings, etc.
bink <- binw/1000
nbp <- adj*2+1

datadir <- paste0(analysisdir,
	"/output/", nbp, "bp_", bink, "k_singletons_", data)

summfile <- paste0(analysisdir, "/summaries/", mac, ".", data, ".summary")
singfile <- paste0(analysisdir, "/singletons/full.singletons")
bindir <- paste0(analysisdir, "/motif_counts/", nbp, "-mers/full")

# Read and preprocess data
full_data <- getData(summfile, singfile, bindir)
gc()

################################################################################
# Drop singletons in individuals with abnormal mutation signatures
################################################################################
i <- 1
cbp <- adj + 1 # I'm guessing cbp = "center base pair"?

bins1Mb <- get_bins(paste0(analysisdir, "/motif_counts/3-mers/full"), "1000kb_full.txt")
bins1Mb$CHR <- gsub("chr", "", bins1Mb$CHR)
# RESULT: CHR, Start, End, Bin, Motif, nMotifs -> how many of motif do we see in each 1000kb window?

sub_mer <- function(motif, n){
	lMotif <- (nchar(motif) - 2) / 2
	center <- ceiling(lMotif / 2)
	SEQA <- substr(motif, center - n, center + n)
	SEQB <- substr(motif, (center * 3) - n, (center * 3) + n)
	return(paste0(SEQA, "(", SEQB, ")"))
}

get_mct_b <- function(bins){
	out <- bins %>%
		group_by(CHR, BIN, Motif) %>%
		summarise(nMotifs=sum(nMotifs))
	return(out)
}

# Count singletons per 3-mer subtype per individual
ind_counts <- full_data$sites %>%
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

# get ids in ped file
ped <- read.table(pedfile, header=F, stringsAsFactors=F)
names(ped)[1] <- "ID"

pheno <- read.table(phenofile, header=T, stringsAsFactors=F)
pheno <- pheno %>%
	dplyr::select(ID=SEQID, Sex, Case.Control, Study, DROP) %>%
	filter(!(is.na(ID)))
# RESULT : ID, Sex, Case.Control, Study, DROP

qplotdat <- read.table(qplotfile, header=T, stringsAsFactors=F)
names(qplotdat)[1] <- "ID"

contam <- read.table(contamfile, header=T, stringsAsFactors=F)
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

# get per-person rates and filter IDs
ind_counts <- merge(ind_counts, mct3, by="Motif") %>%
	filter(ID %in% subjectData$ID) %>%
	mutate(ERV_rel_rate=n/nMotifs,
		subtype=paste0(Type, "_", Motif))

###############################################################################
# get signature loadings
###############################################################################
get_loads <- function(widedat, nmfdat){
	sigloads <- data.frame(subtype=names(widedat),
		sig1=coef(nmfdat)[1,],
		sig2=coef(nmfdat)[2,],
		sig3=coef(nmfdat)[3,]) %>%
		mutate(sig1=sig1/sum(sig1), sig2=sig2/sum(sig2), sig3=sig3/sum(sig3)) %>%
		gather(subtype, value)

	names(sigloads) <- c("subtype", "sig", "value")
	sigloads <- sigloads %>%
		mutate(Category=substr(subtype, 1, 5),
			Sequence=substr(subtype, 7, 14))

	return(sigloads)
}

###############################################################################
# plot signature loadings
###############################################################################
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

###############################################################################
# plot signature contribution across individuals
###############################################################################
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

###############################################################################
# Prepare data and run NMF
###############################################################################
ind_wide <- ind_counts	%>%
	dplyr::select(ID, subtype, ERV_rel_rate) %>%
	spread(subtype, ERV_rel_rate)
ind_wide[is.na(ind_wide)] <- 0
# RESULT: subject x rel_rates_per_3mer

ind_nmf <- nmf(as.matrix(ind_wide[,-c(1)]), 3, nrun=10)
ind_pred <- predict(ind_nmf, "rows", prob=T)

sigloads <- get_loads(ind_wide[,-c(1)], ind_nmf)
plot_loads(sigloads)
ggsave(paste0(analysisdir, "/images/sigloads.png"), width=12, height=6)

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

	plot_ind_sigs(ind_nmf_long)
	ggsave(paste0(analysisdir, "/images/by_ind_all.png"), width=8, height=6)

	plot_ind_sigs2(ind_nmf_long)
	ggsave(paste0(analysisdir, "/images/by_ind_all_by_sig.png"), width=8, height=6)

	###############################################################################
# IDs to keep
###############################################################################
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

ind_nmf_long %>%
  filter(ID %in% keep_ids$ID) %>%
plot_ind_sigs()
ggsave(paste0(analysisdir, "/images/by_ind_post_filter.png"), width=8, height=6)

###############################################################################
# Repeat with kept IDs
###############################################################################
ind_wide2 <- ind_wide %>%
  filter(ID %in% keep_ids$ID)
ind_nmf2 <- nmf(ind_wide2[,-c(1)], 3, nrun=10)
ind_pred2 <- predict(ind_nmf2, "rows", prob=T)

sigloads <- get_loads(ind_wide2[,-c(1)], ind_nmf2)
plot_loads(sigloads)
ggsave(paste0(analysisdir, "/images/sigloads2.png"), width=12, height=6)

nmfdat2 <- data.frame(ID=ind_wide2$ID, basis(ind_nmf2), sig=ind_pred2[[1]])

nmfdat2 <- merge(nmfdat2, subjectData, by="ID") %>%
  mutate(sum=X1+X2+X3,
    sig1=X1/sum,
    sig2=X2/sum,
    sig3=X3/sum)

nmfdat2$ID <- factor(nmfdat2$ID, levels = unique(nmfdat2$ID))

ind_nmf_long <- nmfdat2 %>%
  gather(Signature, prob, sig1:sig3) %>%
  arrange(Study, sum)

plot_ind_sigs(ind_nmf_long)
ggsave(paste0(analysisdir, "/images/by_ind_keep.png"), width=12, height=8)

# test for significant differences between groups
cbind(ind_wide2, g2=nmfdat2$sig) %>%
  mutate(g2=ifelse(g2==1, TRUE, FALSE)) %>%
  mutate(key=paste0(ID, "+", g2)) %>%
  dplyr::select(-c(ID, g2)) %>%
  gather(key2, value, -key) %>% setNames(c("IDf", "key2", "val")) %>%
  tidyr::separate(IDf, into=c("ID", "g2"), sep="[+]") %>%
  mutate(Category=substr(key2, 1, 5),
    Sequence=substr(key2, 7, 14)) %>%
  group_by(Category, Sequence) %>%
  do(tidy(t.test(val~g2, data=.))) %>%
  dplyr::select(Category, Sequence, estimate1, estimate2, p.value) %>%
  filter(p.value<0.05/96)

ggplot()+
  geom_point(data=nmfdat2[nmfdat2$sig==3,],
    aes(x=PC4, y=PC2), colour="gray30", alpha=0.3)+
  geom_point(data=nmfdat2[nmfdat2$sig==1,],
    aes(x=PC4, y=PC2), colour="blue", alpha=0.3)+
  geom_point(data=nmfdat2[nmfdat2$sig==2,],
    aes(x=PC4, y=PC2), colour="yellow", alpha=0.3)+
  theme_bw()
ggsave(paste0(analysisdir, "/images/sig_snp_pcs.png"), width=8, height=8)

full_data$sites <- full_data$sites %>%
	filter(ID %in% keep_ids$ID)

##############################################################################
# Get relative mutation rates per subtype; plot as heatmap
##############################################################################
full_data$aggseq <- get_aggseq(full_data$sites, full_data$mct)
