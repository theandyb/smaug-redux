library(tidyverse)
library(smaug)
library(yaml)

args <- yaml.load_file("/net/snowwhite/home/beckandy/research/smaug-redux/_config.yaml")
attach(args)

bink <- binw/1000
nbp <- adj*2+1

datadir <- paste0(analysisdir, "/output/", nbp, "bp_", bink, "k_singletons_", data)

summfile <- paste0(analysisdir, "/summaries/", mac, ".", data, ".summary")
singfile <- paste0(analysisdir, "/singletons/full.singletons")
maskfile <- paste0(analysisdir, "/reference_data/testmask2.bed")

getData2 <- function(summfile, singfile, maskfile, binw){
	cat("Reading summary file:", summfile, "...\n")
	sites <- read_tsv(summfile)
	sites$Category2 <- sites$Category
	sites$Category <- gsub("cpg_", "", sites$Category)
	sites$BIN <- ceiling(sites$POS / binw)
	sites$MASK <- binaryCol(sites, maskfile)

	cat("Annotating with sample ID...\n")
	inds <- read_tsv(singfile)
	names(inds) <- c("CHR", "POS", "S", "ALT", "ID")
	inds <- inds %>% dplyr::select(CHR, POS, ALT, ID)
	sites <- merge(sites, inds, by=c("CHR", "POS", "ALT"))

	return(sites)
}

df <- getData2(summfile, singfile, maskfile, binw)
summfile2 <- paste0(summfile, ".annotated")
print("Data loaded, writing output")
write.table(df, summfile2, quote = FALSE, row.names = FALSE, sep = "\t")
