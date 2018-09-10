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
  print(str_c("Processing file: ", pred.name))
  # load the data
  dat <- read_tsv(pred.name, col_names = FALSE)
  chr <- pred.name %>% str_extract("chr[0-9]+") %>% str_extract("[0-9]+") %>% as.numeric()
  type <- str_extract(pred.name, "[A,C,G,T]{2}_[A,C,G,T]{2}")
  ret <- dat %>% 
          mutate(chr = chr,
                 type = type,
                 min.rate = min(X3),
                 max.rate = max(X3),
                 median.rate = median(X3),
                 var.rate = var(X3),
                 nSites = n(),
                 nMutableSites = sum(X3 > 0.5)
                )
  return(ret)
}

res <- mclapply(pred.files, get_sum, mc.cores = detectCores())

res.table <- bind_rows(res)

write.csv(res.table, file = "explore.csv", row.names=FALSE)
