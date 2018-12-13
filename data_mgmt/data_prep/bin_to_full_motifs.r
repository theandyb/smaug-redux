library(tidyverse)
library(here)
library(yaml)

args <- yaml.load_file(here("_config.yaml"))
results_dir <- paste0(args$analysisdir, "/motif_counts/3-mers/full")

chromosomes <- seq(1,22)

for(chr in chromosomes){
    binFile <- paste0(results_dir, "/", "chr", chr, ".3-mer_motifs_1000kb_full.txt")
    df <- read_tsv(binFile)
    df2 <- df %>% group_by(Motif, CHR) %>% summarize(nMotifs = sum(nMotifs)) %>% select(CHR, Motif, nMotifs)
    binFile <- paste0(results_dir, "/", "chr", chr, ".3-mer_motifs_full.txt")
    write.table(df2, sep="\t", outFile, quote = FALSE, row.name = FALSE)
}