library("tidyverse")
library("yaml")
initial.options <- commandArgs(trailingOnly = FALSE)
script.basename <- sub("--file=", "", initial.options[grep("--file=", initial.options)]) %>%
  dirname
setwd(script.basename)
args <- yaml.load_file("../_config.yaml")
print(args)
