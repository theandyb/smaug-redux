# What this is? Everything I've run in order to run the validation code
# How Jed wrote this: run mod_shell.R again after fitting logistic regression
# I.e. a lot of redundancy

# Goal: ID what needed to be saved off in first run that wasn't

###################################################################################
# Step 1: Read and pre-process data
###################################################################################

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

##########################################################################################
# Step 2: Remove subjects with abnormal mutation profiles (NMF)
##########################################################################################

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
