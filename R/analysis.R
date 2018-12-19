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
packages <- c("broom", "RColorBrewer", "MASS", "boot", "speedglm",
	"psych", "lmtest", "fmsb", "hexbin", "cowplot", "grid", "gtable", "gridExtra",
	"openxlsx", "Biostrings", "svglite", "NMF", "emdbook", "tidyverse")
invisible(sapply(packages, function(x)
	suppressMessages(require(x, character.only = TRUE))))

# Define additional variables for cleaner strings, etc.
bink <- binw/1000
nbp <- adj*2+1

datadir <- paste0(analysisdir,
	"/output/", nbp, "bp_", bink, "k_singletons_", data)

summfile <- paste0(analysisdir, "/summaries/", mac, ".", data, ".summary")
singfile <- paste0(analysisdir, "/singletons/full.singletons")
bindir <- paste0(analysisdir, "/motif_counts/", nbp, "-mers/full")
maskfile <- paste0(analysisdir, "/reference_data/testmask2.bed")

# Read and preprocess data
full_data <- getData(summfile, singfile, bindir, maskfile, binw)
gc()

################################################################################
# Drop singletons in individuals with abnormal mutation signatures
################################################################################
source('R/subject_filter.r')
binDir <- paste0(analysisdir, "/motif_counts/3-mers/full")
binPattern <- "1000kb_full.txt"

subjectData <- get_ind_data(pedfile, phenofile, qplotfile, contamfile, vcfastfile, pcfile)
ind_wide <- get_nmf_input(sites, binDir, binPattern, subjectData)
keep_list <- nmf_analysis(ind_wide, subjectData, analysisdir, outPlots = FALSE)


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
cbp <- 5

ratelist <- list()
testlist <- list()
modlist <- list()
for(j in 2:5){
	i <- j-1

	# using nMotifs in the full file results in incorrect counts when collapsing
	# shorter motifs. Here we merge with binned outputs obtained from running
	# the data_pipeline/augment_summary.pl script with shorter motifs to ensure
	# proper counting of mutable sites
	nbptmp <- i*2+1
	bindir <- paste0(analysisdir, "/motif_counts/", nbptmp, "-mers/full")
  p1 <- "motifs_full.txt"
	gpdat <- full_data$aggseq %>%
		mutate(Type=gsub("cpg_", "", Category2),
			SEQA=substr(Motif, cbp-i, cbp+i),
			SEQB=substr(Motif, cbp*3-i, cbp*3+i),
			Motif=paste0(SEQA, "(", SEQB, ")"))

	if(i==0){
		gpdat <- gpdat %>%
			dplyr::select(Type, Motif, nERVs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs))

		mcfile <- paste0(analysisdir, "/output/", nbptmp, "bp_final_rates.txt")
		mcount <- read.table(mcfile, header=T, stringsAsFactors=F)
		mcount <- mcount %>%
			mutate(Motif=ifelse(grepl("^A", Type),
				"A(T)",
				"C(G)")) %>%
			dplyr::select(Type, Motif, nMotifs)

		gpdat <- merge(gpdat, mcount, by=c("Type", "Motif")) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
	}
	else if(i>0 & i<4){
 		gpdat<- gpdat %>%
			dplyr::select(Type, Motif, nERVs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs))

		bins <- get_bins(bindir, p1)
		mcount <- get_mct(bins)

		gpdat <- merge(gpdat, mcount, by=c("Motif")) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)

		cat("Plotting heatmap panels...\n")
		for(k in 1:6){
			categ <- orderedcats[k]

			p1 <- rrheat2(gpdat[gpdat$Type==categ,], i)
			p1a <- p1+theme(legend.position="none")

			setEPS()
			postscript(paste0(analysisdir, "/images/", categ, "_", nbptmp, "bp_heatmap.eps"),
			height=5, width=5)
			pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
			grid.draw(ggplotGrob(p1a))
			dev.off()
		}

		cat("Trimming panel whitespace...\n")
		# trim whitespace on panels with imagemagick mogrify
		trimcmd <- paste0("mogrify -trim ",
			analysisdir, "/images/*", nbptmp, "bp_heatmap.eps")
		system(trimcmd)

		# extract legend
		legend <- get_legend(p1)
		setEPS()
		postscript(paste0(analysisdir, "/images/heatmap_legend.eps"),
			height=8, width=3)
		grid.draw(legend)
		dev.off()

	}
	else { # don't draw heatmap for 9-mers; use original motif counts
		gpdat <- gpdat %>%
			dplyr::select(Type, Motif, nERVs, nMotifs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs), nMotifs=sum(nMotifs)) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
	}

	ratelist[[j]] <- gpdat

	# Test for heterogeneity among subtypes sharing same (K-2)-mer parent
	if(i>0){
		parentdat <- gpdat %>%
			# dplyr::select(Type, Motif, nERVs, nMotifs, rel_prop) %>%
			# filter(nERVs > 20) %>%
			mutate(SEQA=substr(Motif, 2, nbptmp-1),
				SEQB=substr(Motif, nbptmp+3, nchar(Motif)-2),
				MotifP=paste0(SEQA, "(", SEQB, ")")) %>%
			group_by(Type, MotifP) %>%
			arrange(Type, MotifP) %>%
			mutate(exp=sum(nERVs)/sum(nMotifs)*nMotifs,
				p=exp/sum(exp),
				n=n(),
				b1=substr(Motif,1,1),
				b2=substr(Motif,nbptmp,nbptmp)) %>%
			filter(n==16)

		moddat <- parentdat %>%
			do(tidy(glance(lm(ERV_rel_rate ~ b1+b2, data=.))))
		modlist[[i]] <- moddat

		hettests <- parentdat %>%
			summarise(pval=chisq.test(nERVs, p=p)$p.value) %>%
			ungroup() %>%
			mutate(fdr=p.adjust(pval, method="fdr"))
		testlist[[i]] <- hettests
	}

	write.table(gpdat[,c("Type", "Motif", "nERVs", "nMotifs", "ERV_rel_rate")],
		paste0(analysisdir, "/output/rates/", nbptmp, "bp_final_rates2.txt"),
		col.names=T, row.names=F, quote=F, sep="\t")
}

##############################################################################
# Compare full data to strict ancestral data
##############################################################################
i <- 3
nbp2 <- 7

# correct (HC)
sites_c_hc <- full_data$sites %>%
  filter(!(tolower(ALT)==tolower(AA)))

sites_c_hc <- sites_c_hc %>%
	mutate(SEQA=substr(Motif, cbp-i, cbp+i),
		SEQB=substr(Motif, cbp*3-i, cbp*3+i),
		Motif=paste0(SEQA, "(", SEQB, ")"))

prop <- nrow(sites_c_hc)/nrow(full_data$sites)

bindir <- paste0(analysisdir, "/motif_counts/", nbp2, "-mers/full")
p1 <- "motifs_full.txt"
bins <- get_bins(bindir, p1)
mct <- get_mct(bins)
ancaggseq <- get_aggseq(sites_c_hc, mct)

maskdir <- paste0(analysisdir, "/motif_counts/", nbp2, "-mers/mask")
p1 <- "motifs_mask.txt"
maskbins <- get_bins(maskdir, p1)
maskmct <- get_mct(maskbins)

ancgpdat <- ancaggseq %>%
  mutate(Type=gsub("cpg_", "", Category2)) %>%
  dplyr::select(Type, Motif,
		nERVs_anc=nERVs, nMotifs_anc=nMotifs,
		ERV_rel_rate_anc=rel_prop)

ancgpdat <- merge(ratelist[[4]], ancgpdat, by=c("Type", "Motif"))

test_anc <- ancgpdat %>%
  group_by(Type, Motif) %>%
  do(tidy(prop.test(c(.$nERVs, .$nERVs_anc), c(.$nMotifs, .$nMotifs_anc*prop)))) %>%
	dplyr::select(Type, Motif, estimate1, estimate2, statistic, p.value)

test_anc <- merge(test_anc, ancgpdat, by=c("Type", "Motif"))

sig_anc <- test_anc %>%
  filter(p.value<0.05/nrow(test_anc)) %>%
  mutate(prop=estimate1/estimate2,
    gt10=ifelse(abs(log(prop))>log(1.2), TRUE, FALSE)) %>%
	filter(gt10==TRUE) %>%
	arrange(desc(prop))


if(build_logit){
	cat("Preparing data for logistic regression model...\n")

	i<-3
	for(chr in 1:22){
		posfile <- paste0(analysisdir,
		"/output/logmod_data/chr", chr, "_sites.txt")
		dat <- full_data$sites %>%
			# summfile1 <- sites %>%
			# filter(Category==categ) %>%
			filter(CHR==chr) %>%
			mutate(Type=gsub("cpg_", "", Category2),
				SEQA=substr(Motif, cbp-i, cbp+i),
				SEQB=substr(Motif, cbp*3-i, cbp*3+i),
				Sequence=paste0(SEQA, "(", SEQB, ")")) %>%
			dplyr::select(CHR, POS, Sequence, Type) %>%
			mutate(mut=1) %>%
			spread(Type, mut, fill=0)

		write.table(dat, posfile, col.names=F, row.names=F, quote=F, sep="\t")
	}
}
