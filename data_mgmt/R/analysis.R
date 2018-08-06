#!/usr/bin/RScript

cat("Loading packages...")
install.packages("yaml", quiet=TRUE)
require(yaml)

args <- yaml.load_file("./_config.yaml")
attach(args)

# add user libpath, if it exists (useful if running on a different node)
.libPaths(c(.libPaths(), args$libpath))

require(devtools)
install_github('carjed/smaug', quiet=TRUE)
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

datadir <- paste0(parentdir,
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
