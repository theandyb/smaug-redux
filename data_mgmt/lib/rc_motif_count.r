library(tidyverse)
library(stringi)
chromosomes <- seq(1,22)
results_dir <- "/net/snowwhite/home/beckandy/research/smaug-redux/motif_counts/3-mers/full"

for(chr in chromosomes){
    cntFile <- paste0(results_dir, "/", "chr", chr, "_bins.tsv")
    df <- read_tsv(cntFile, col_names = c('Motif', 'nMotif', 'BIN'))
    df <- df %>% mutate(RC = stri_reverse(chartr("ACGT","TGCA", Motif)))
    df <- inner_join(select(df, Motif, RC, BIN, nMotif), 
        select(df, RC, BIN, nMotif),
        by = c("Motif" = "RC", "BIN" = "BIN"))
    df['nMotifs'] <- df$nMotif.x + df$nMotif.y
    df <- df %>% mutate(Motif = paste0(Motif, "(", RC ,")")) %>%
        select(Motif, BIN, nMotifs)
    df['CHR'] <- chr
    outFile <- paste0(results_dir, "/", "chr", chr, ".3-mer_motifs_1000kb_full.txt")
    write.table(df, outFile, quote = FALSE, row.name = FALSE)
}
