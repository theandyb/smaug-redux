library(tidyverse)
library(stringi)
library(here)
library(yaml)

args <- yaml.load_file(here("_config.yaml"))
attach(args)
k <- "7"
chromosomes <- seq(1,22)

results_dir <- paste0(analysisdir, "/motif_counts/",k,"-mers/full")

for(chr in chromosomes){
    cntFile <- paste0(results_dir, "/", "chr", chr, "_bins.tsv")
    df <- read_tsv(cntFile, col_names = c('Motif', 'nMotifs', 'BIN'))
    df <- df %>% mutate(RC = stri_reverse(chartr("ACGT","TGCA", Motif)))
    df <- inner_join(select(df, Motif, RC, BIN, nMotifs), 
        select(df, RC, BIN, nMotifs),
        by = c("Motif" = "RC", "BIN" = "BIN"))
    df['nMotifs'] <- df$nMotifs.x + df$nMotifs.y
    df <- df %>% mutate(Motif = paste0(Motif, "(", RC ,")")) %>%
        select(Motif, BIN, nMotifs)
    df['CHR'] <- chr
    outFile <- paste0(results_dir, "/", "chr", chr, ".",k,"-mer_motifs_1000kb_full.txt")
    write.table(df, outFile, sep="\t", quote = FALSE, row.name = FALSE)
}
