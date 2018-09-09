library("tidyverse")
library("stringr")
library("yaml")
library("parallel")
initial.options <- commandArgs(trailingOnly = FALSE)
script.basename <- sub("--file=", "", initial.options[grep("--file=", initial.options)]) %>%
  dirname
setwd(script.basename)

args <- yaml.load_file("../_config.yaml")
attach(args)

pred.files.dir <- paste0(analysisdir, "/output/predicted")
setwd(pred.files.dir)

pred.files <- Sys.glob("chr*.txt")
#for(pred.name in pred.files){
get_sum <- function(pred.name){
  # What do we want to get? Min, Max, Var, Number > 0.5
  # load the data
  dat <- read_tsv(pred.name, col_names = FALSE)
  chr <- pred.name %>% str_extract("chr[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric
  type <- str_extract(pred.name, "[A,C,G,T]{2}_[A,C,G,T]{2}")
  return tibble(chr = chr, type = type, min = min(dat$X3), max = max(dat$X3), var = var(dat$X3), highly = sum(dat$X3 > 0.5))
}

res <- mclapply(pred.files, get_sum, mc.cores = detectCores())

res.table <- bind_rows(res)

write.csv(res.table, file = "explore.csv")
