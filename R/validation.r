# Libraries and global variables
library(smaug)
library(yaml)
library(lmtest)
library(tidyverse)
library(magrittr)

args <- yaml.load_file("./_config.yaml")
attach(args)

# Read and process the data
avrates <- read_excel(paste0(parentdir, "/reference_data/AV_rates.xlsx"), sheet=10)
validation_file <- paste0(parentdir, "/output/predicted/validation_sites.txt")
input_sites <- read_tsv(validation_file, col_names = c("CHR", "POS", "MU", "OBS", "Category", "SEQ", "ID"))

input_sites %<>% filter(MU > 0)
input_dnms <- input_sites %>% filter(ID != "all") %>% mutate(SIM = "a", SIMOBS = 0)
input_sites %<>% filter(ID == "all")

m500k <- validationPipe(input_sites, input_dnms, 500000, ratelist, avrates)
m1m <- validationPipe(input_sites, input_dnms, 1000000, ratelist, avrates)
m2m <- validationPipe(input_sites, input_dnms, 2000000, ratelist, avrates)
m3m <- validationPipe(input_sites, input_dnms, 3000000, ratelist, avrates)