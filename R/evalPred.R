library("tidyverse")
library("yaml")
initial.options <- commandArgs(trailingOnly = FALSE)
script.basename <- sub("--file=", "", initial.options[grep("--file=", initial.options)]) %>%
  dirname
setwd(script.basename)

args <- yaml.load_file("../_config.yaml")
attach(args)

pred.files.dir <- paste0(analysisdir, "/output/predicted")
setwd(pred.files.dir)

pred.files <- Sys.glob("chr*.txt")
print(pred.files)
