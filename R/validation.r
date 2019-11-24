# Libraries and global variables
library(smaug)
library(yaml)
library(lmtest)
library(smaug)
library(tidyverse)
library(broom)
library(readxl)
library(magrittr)

args <- yaml.load_file("./_config.yaml")
attach(args)

# We need the ratelist from the kmer analysis
sitesFile <- paste0(analysisdir, "/output/intermediates/filtered.sites.csv.gz")

bink <- binw/1000
nbp <- adj*2+1

bindir <- paste0(analysisdir, "/motif_counts/", nbp, "-mers/full")

full_data <- list()
full_data$sites <- read_csv(sitesFile)
full_data$bins <- get_bins(bindir, "motifs_full.txt")
full_data$mct <- get_mct(full_data$bins)

full_data$aggseq <- get_aggseq(full_data$sites, full_data$mct)
cbp <- 5
p1 <- "motifs_full.txt"

# Generate r5m (from AV_rates.r)
avrates <- read_excel(paste0(analysisdir, "/reference_data/AV_rates.xlsx"), sheet=10)
names(avrates) <- c("Motif", "alt", "afr", "asn", "eur", "refrev", "altrev")
rates7 <- ratelist[[4]]

avrates$CAT <- paste0(substr(avrates$Motif,4,4), substr(avrates$alt, 4,4))
avrates$Category[avrates$CAT=="AC" | avrates$CAT=="TG"] <- "AT_CG"
avrates$Category[avrates$CAT=="AG" | avrates$CAT=="TC"] <- "AT_GC"
avrates$Category[avrates$CAT=="AT" | avrates$CAT=="TA"] <- "AT_TA"
avrates$Category[avrates$CAT=="GA" | avrates$CAT=="CT"] <- "GC_AT"
avrates$Category[avrates$CAT=="GC" | avrates$CAT=="CG"] <- "GC_CG"
avrates$Category[avrates$CAT=="GT" | avrates$CAT=="CA"] <- "GC_TA"

avrates$Motif <- paste0(avrates$Motif, "(", avrates$refrev, ")")

avrates <- avrates %>%
  dplyr::select(Type=Category, Motif, eur)

r5m <- merge(rates7, avrates, by=c("Type", "Motif"))

r5m$Category2 <- ifelse(substr(r5m$Motif,4,5)=="CG",
  paste0("cpg_",r5m$Type),
  r5m$Type)

r5m <- r5m %>%
  mutate(Category2 = plyr::mapvalues(Category2, orderedcats1, orderedcats2))

r5m$Category2 <- factor(r5m$Category2, levels=orderedcats2)

r5m$prop_diff <- (r5m$eur/(mean(r5m$eur)/mean(r5m$ERV_rel_rate)))/r5m$ERV_rel_rate
r5m$prop_diff4 <- r5m$prop_diff
r5m$prop_diff4[r5m$prop_diff< 0.5] <- 0.5
r5m$prop_diff4[r5m$prop_diff>2] <- 2

r5m$v2 <- substr(r5m$Motif,1,3)
r5m$v2a <- as.character(lapply(as.vector(r5m$v2), reverse_chars))
r5m$v2a <- factor(r5m$v2a)
r5m$v3 <- substr(r5m$Motif,3+2,3*2+1)
rm(rates7)

kmerAnalysis <- function(aggseq, analysisdir, p1, cbp, i, plotOutput = TRUE) {
	testlist <- list()
	modlist <- list()
	j <- i+1
	nbptmp <- i*2+1
	bindir <- paste0(analysisdir, "/motif_counts/", nbptmp, "-mers/full")
    
	gpdat <- aggseq %>%
		mutate(Type=gsub("cpg_", "", Category2),
			SEQA=substr(Motif, cbp-i, cbp+i),
			SEQB=substr(Motif, cbp*3-i, cbp*3+i),
			Motif=paste0(SEQA, "(", SEQB, ")"))

	if(i>0 & i<4){
 		gpdat<- gpdat %>%
			dplyr::select(Type, Motif, nERVs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs))

		bins <- get_bins(bindir, p1)
		mcount <- get_mct(bins)

		gpdat <- merge(gpdat, mcount, by=c("Motif")) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
		if (plotOutput){
			cat("Plotting heatmap panels...\n")
			for(k in 1:6){
				categ <- orderedcats[k]

				p1 <- rrheat2(gpdat[gpdat$Type==categ,], i)
				p1a <- p1+theme(legend.position="none")

				setEPS()
				postscript(paste0(analysisdir, "/images/", categ, "_", nbptmp, "bp_heatmap.eps"),
				height=5, width=5)
				pushViewport(viewport(width=unit(5, "in"), height=unit(5, "in")))
				grid.draw(ggplotGrob(p1a))
				dev.off()
			}

			cat("Trimming panel whitespace...\n")
			# trim whitespace on panels with imagemagick mogrify
			trimcmd <- paste0("mogrify -trim ",
				analysisdir, "/images/*", nbptmp, "bp_heatmap.eps")
			system(trimcmd)

			# extract legend
			legend <- get_legend(p1)
			setEPS()
			postscript(paste0(analysisdir, "/images/heatmap_legend.eps"),
				height=8, width=3)
			grid.draw(legend)
			dev.off()
		}

	}
	else { # don't draw heatmap for 9-mers; use original motif counts
		gpdat <- gpdat %>%
			dplyr::select(Type, Motif, nERVs, nMotifs) %>%
			group_by(Type, Motif) %>%
			summarise(nERVs=sum(nERVs), nMotifs=sum(nMotifs)) %>%
			mutate(ERV_rel_rate=nERVs/nMotifs)
	}

	ratelist <- gpdat

	# Test for heterogeneity among subtypes sharing same (K-2)-mer parent
	if(i>0){
		parentdat <- gpdat %>%
			mutate(SEQA=substr(Motif, 2, nbptmp-1),
				SEQB=substr(Motif, nbptmp+3, nchar(Motif)-2),
				MotifP=paste0(SEQA, "(", SEQB, ")")) %>%
			group_by(Type, MotifP) %>%
			arrange(Type, MotifP) %>%
			mutate(exp=sum(nERVs)/sum(nMotifs)*nMotifs,
				p=exp/sum(exp),
				n=n(),
				b1=substr(Motif,1,1),
				b2=substr(Motif,nbptmp,nbptmp)) %>%
			filter(n==16)

		moddat <- parentdat %>%
			do(broom::glance(lm(ERV_rel_rate ~ b1+b2, data=.)))
		modlist <- moddat

		hettests <- parentdat %>%
			summarise(pval=chisq.test(nERVs, p=p)$p.value) %>%
			ungroup() %>%
			mutate(fdr=p.adjust(pval, method="fdr"))
		testlist <- hettests
	}

	write.table(gpdat[,c("Type", "Motif", "nERVs", "nMotifs", "ERV_rel_rate")],
		paste0(analysisdir, "/output/rates/", nbptmp, "bp_final_rates2.txt"),
		col.names=T, row.names=F, quote=F, sep="\t")
	return(list(modlist, ratelist, testlist))
}

ratelist <- list()
testlist <- list()
modlist <- list()

for(i in 1:4){
	res <- kmerAnalysis(full_data$aggseq, analysisdir, p1, cbp, i, plotOutput = FALSE)
	ratelist[[i+1]] <- res[[2]] 
	testlist[[i]] <- res[[3]]
	modlist[[i]] <- res[[1]]
	rm(res)
}

res <- kmerAnalysis(full_data$aggseq, analysisdir, p1, cbp, 0, plotOutput = FALSE)
ratelist[[1]] <- res[[2]] 
rm(res)

# Read and process the data
validation_file <- paste0(analysisdir, "/output/predicted/validation_sites.txt")
input_sites <- read_tsv(validation_file, col_names = c("CHR", "POS", "MU", "OBS", "Category", "SEQ", "ID"))

input_sites %<>% filter(MU > 0)
input_dnms <- input_sites %>% filter(ID != "all") %>% mutate(SIM = "a", SIMOBS = 0)
input_sites %<>% filter(ID == "all")

m500k <- validationPipe(input_sites, input_dnms, 500000, ratelist, avrates, rm5)
m1m <- validationPipe(input_sites, input_dnms, 1000000, ratelist, avrates, rm5)
m2m <- validationPipe(input_sites, input_dnms, 2000000, ratelist, avrate, rm5s)
m3m <- validationPipe(input_sites, input_dnms, 3000000, ratelist, avrates, rm5)