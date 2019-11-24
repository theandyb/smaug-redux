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

orderedcats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA", "cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA")
orderedcats1 <- c("AT_GC", "GC_AT", "cpg_GC_AT", "AT_CG", "GC_CG", "cpg_GC_CG", "AT_TA", "GC_TA", "cpg_GC_TA")
orderedcats2 <- c("A>G", "C>T", "CpG>TpG", "A>C", "C>G", "CpG>GpG", "A>T", "C>A", "CpG>ApG")

datadir <- paste0(analysisdir,
	"/output/", nbp, "bp_", bink, "k_singletons_", data)

summfile <- paste0(analysisdir, "/summaries/", mac, ".", data, ".summary")
singfile <- paste0(analysisdir, "/singletons/full.singletons")
bindir <- paste0(analysisdir, "/motif_counts/", nbp, "-mers/full")
maskfile <- paste0(analysisdir, "/reference_data/testmask2.bed")

sitesFile <- paste0(analysisdir, "/output/intermediates/filtered.sites.csv.gz")

if (!file.exists(sitesFile)){
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
	ind_wide <- get_nmf_input(full_data$sites, binDir, binPattern, subjectData)
	keep_ids <- nmf_analysis(ind_wide, subjectData, analysisdir, outPlots = FALSE)

	if (!dir.exists(paste0(analysisdir, "/output/intermediates"))){
		dir.create(paste0(analysisdir, "/output/intermediates"))
	}

	write.csv(keep_ids, file = paste0(analysisdir, "/output/intermediates/keep_ids.csv"),
		quote = FALSE, row.names = FALSE)

	full_data$sites <- full_data$sites %>%
		filter(ID %in% keep_ids$ID)

	write.csv(full_data$sites, file=gzfile(sitesFile), quote = FALSE, row.names = FALSE)
} else {
	full_data <- list()
	full_data$sites <- read_csv(sitesFile)
	full_data$bins <- get_bins(bindir, "motifs_full.txt")
	full_data$mct <- get_mct(full_data$bins)

	keep_ids <- read_csv(paste0(analysisdir, "/output/intermediates/keep_ids.csv"))
}
##############################################################################
# Get relative mutation rates per subtype; plot as heatmap
##############################################################################
full_data$aggseq <- get_aggseq(full_data$sites, full_data$mct)
cbp <- 5
p1 <- "motifs_full.txt"

kmerAnalysis <- function(aggseq, analysisdir, p1, cbp, i, plotOutput = TRUE) {
	testlist <- list()
	modlist <- list()
	j <- i+1
	nbptmp <- i*2+1
	bindir <- paste0(analysisdir, "/motif_counts/", nbptmp, "-mers/full")
    
	gpdat <- aggseq %>%
		mutate(Type=gsub("cpg_", "", Category2),
			SEQA=substr(Motif, cbp-i, cbp+i),
			SEQB=substr(Motif, cbp*3-i, cbp*3+i),
			Motif=paste0(SEQA, "(", SEQB, ")"))

	if(i>0 & i<4){
 		gpdat<- gpdat %>%
			dplyr::select(Type, Motif, nERVs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs))

		bins <- get_bins(bindir, p1)
		mcount <- get_mct(bins)

		gpdat <- merge(gpdat, mcount, by=c("Motif")) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
		if (plotOutput){
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
			do(glance(lm(ERV_rel_rate ~ b1+b2, data=.)))
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
	return(list(modlist, ratelist, testlist))
}

ratelist <- list()
testlist <- list()
modlist <- list()

for(i in 1:4){
	res <- kmerAnalysis(full_data$aggseq, analysisdir, p1, cbp, i, ratelist, testlist, modlist, plotOutput = FALSE)
	modlist[[i]] <- res[[1]]
	ratelist[[i+1]] <- res[[2]] 
	testlist[[i]] <- res[[3]]
	rm(res)
}

res <- kmerAnalysis(full_data$aggseq, analysisdir, p1, cbp, 0, ratelist, testlist, modlist, plotOutput = FALSE)
ratelist[[1]] <- res[[2]] 
rm(res)

# Prep data for logistic regression models

logmodData <- function(sites, chr, outFile, cbp, i){
	dat <- full_data$sites %>%
			filter(CHR==chr) %>%
			mutate(Type=gsub("cpg_", "", Category2),
				SEQA=substr(Motif, cbp-i, cbp+i),
				SEQB=substr(Motif, cbp*3-i, cbp*3+i),
				Sequence=paste0(SEQA, "(", SEQB, ")")) %>%
			dplyr::select(CHR, POS, Sequence, Type) %>%
			mutate(mut=1) %>%
			spread(Type, mut, fill=0)
	write_tsv(dat, outFile, col_names = FALSE, quote_escape = FALSE)
}

if(build_logit){
	cat("Preparing data for logistic regression model...\n")

	i<-3
	for(chr in 1:22){
		cat(paste0("Chromosome ", chr, "...\n"))
		posfile <- paste0(analysisdir, "/output/logmod_data/chr", chr, "_sites.txt")
		logmodData(full_data$sites, chr, posfile, cbp, i)
	}
}

# Code for evaluating the rates