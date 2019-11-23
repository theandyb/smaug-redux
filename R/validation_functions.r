#' Validation models
#'
#' Runs a set of logistic regression models for the validation step of the pipeline
#' 
#' @param data
#' @param group set to "COMBINED" to run combined models
#'
#' @return List with glm objects
#'
#' @export
runDNMLogit <- function(data, group){
	models <- list()
	if(group=="COMBINED"){
		# Nested models (1-5)
		models$mod_1mers_n <- glm(OBS~MU_1,
			data=data, family=binomial())
		models$mod_3mers_n <- glm(OBS~MU_1+resid3,
			data=data, family=binomial())
		models$mod_5mers_n <- glm(OBS~MU_1+resid3+resid5,
			data=data, family=binomial())
		models$mod_7mers_n <- glm(OBS~MU_1+resid3+resid5+resid7,
			data=data, family=binomial())
		models$mod_7mers_features_n <- glm(OBS~MU_1+resid3+resid5+resid7+residL,
			data=data, family=binomial())

		# Non-nested models (6-8)
		models$mod_1mers <- glm(OBS~MU_1, data=data, family=binomial())
		models$mod_3mers <- glm(OBS~MU_3, data=data, family=binomial())
		models$mod_5mers <- glm(OBS~MU_5, data=data, family=binomial())

		# Comparison of 7-mer models (9-13)
		models$mod_7mers <- glm(OBS~MU_7, data=data, family=binomial())
		models$mod_7mers_ERVs_down <- glm(OBS~MU_7D, data=data, family=binomial())
		models$mod_7mers_MAC10 <- glm(OBS~MU_7C, data=data, family=binomial())
		models$mod_7mers_AV <- glm(OBS~MU_7A, data=data, family=binomial())
		models$mod_7mers_masked <- glm(OBS~MU_7M, data=data, family=binomial())

		# Non-nested 7-mer+features model (14)
		models$mod_7mers_features <- glm(OBS~MU, data=data, family=binomial())
		models$mod_7mers_adj <- glm(OBS~MU_7adj, data=data, family=binomial())
		models$mod_7mers_anc <- glm(OBS~MU_7AN, data=data, family=binomial())
		# models$mod_7mers_sig <- glm(OBS~MU_7+X1+X2+X3+X4+X5, data=data, family=binomial())
		# models$mod_9mers <- glm(OBS~MU_9, data=data, family=binomial())

		models$mod_int <- glm(OBS~1,
			data=data, family=binomial())
	} else {
		# Nested models (1-4)
		models$mod_3mers_n <- glm(OBS~MU_3,
			data=data, family=binomial())
		models$mod_5mers_n <- glm(OBS~MU_3+resid5,
			data=data, family=binomial())
		models$mod_7mers_n <- glm(OBS~MU_3+resid5+resid7,
			data=data, family=binomial())
		models$mod_7mers_features_n <- glm(OBS~MU_3+resid5+resid7+residL,
			data=data, family=binomial())

		# Non-nested models (5-6)
		models$mod_3mers <- glm(OBS~MU_3, data=data, family=binomial())
		models$mod_5mers <- glm(OBS~MU_5, data=data, family=binomial())

		# Comparison of 7-mer models (7-11)
		models$mod_7mers <- glm(OBS~MU_7, data=data, family=binomial())
		models$mod_7mers_ERVs_down <- glm(OBS~MU_7D, data=data, family=binomial())
		models$mod_7mers_MAC10 <- glm(OBS~MU_7C, data=data, family=binomial())
		models$mod_7mers_AV <- glm(OBS~MU_7A, data=data, family=binomial())
		models$mod_7mers_masked <- glm(OBS~MU_7M, data=data, family=binomial())

		# Non-nested 7-mer+features model (12)
		models$mod_7mers_features <- glm(OBS~MU, data=data, family=binomial())
		models$mod_7mers_anc <- glm(OBS~MU_7AN, data=data, family=binomial())
		# models$mod_7mers_sig <- glm(OBS~MU_7+X1+X2+X3+X4+X5, data=data, family=binomial())
		# models$mod_9mers <- glm(OBS~MU_9, data=data, family=binomial())
	}

	return(models)
}

#' Get Nagelkerke's R^2 From Model List
#'
#' Parse over lis of glm objects and return Nagelkerke's R^2 for each
#'
#' @param models A list object containing glm objects
#' 
#' @return A vector of R^2 values
#'
#' @importFrom fmsb NagelkerkeR2
#' @export
getNR2 <- function(models){
	unlist(lapply(models, function(x)
		fmsb::NagelkerkeR2(x)))[seq(2, 2 * length(models), by = 2)]
}

#' Get AIC values for model objects in a list
#'
#' Parse over lis of glm objects and return AIC for each
#'
#' @param models A list object containing glm objects
#'
#' @return A vector containing AIC values
#'
#' @importFrom stats AIC
#' @export
getAIC <- function(models){
	unlist(lapply(models, function(x)
		AIC(x)))[1:length(models)]
}

#' Simulate mutation data
#' 
#' This function returns simulated mutation data 
#'
#' @param data
#' @param EST
#' @param nobs
#' @param chunksize
#' 
#' @return data.frame containing simulated mutation data
#' @export
simMu <- function(data, EST, nobs, chunksize=50000){
	mutated <- data.frame()
	N <- nrow(data)
	while(nrow(mutated) < nobs){
		rowind <- sample(N, chunksize)
		batch <- data[rowind,]
		mu <- batch[,c(EST)]

		batch$SIMOBS <- sapply(mu, function(x) rbinom(1,1,x))
		mutated <- rbind(mutated, batch[batch$SIMOBS == 1,])
	}

	mutated <- mutated[sample(nrow(mutated), nobs),] %>%
		mutate(OBS=0, SIM="b")
	return(mutated)
}

#' Sample non-mutated sites
#'
#' Samples a given number of non-mutated sites from a data set
#' 
#' @param data
#' @param nsites
#'
#' @return data.frame with a subsample of the input data set
#' @import magrittr
#' @export
buildValidationData <- function(data, nsites){
    set.seed(nsites)
    outdat <- data[sample(nrow(data), nsites),] %>%
        mutate(SIM="ab", SIMOBS=0)

    return(outdat)
}

#' Merge rate estimates with validation data
#' 
#' A longer description once you've figured out how this works
#' @import magrittr
#' @export
mergeRates <- function(chrp_c, avrates, ratelist){
    adj2 <- 3

	rates7 <- ratelist[[4]] %>%
		mutate(SEQ7=substr(Motif, 1, 7)) %>%
		dplyr::select(Category=Type, SEQ7, MU_7=ERV_rel_rate)

	rates5 <- ratelist[[3]] %>%
		mutate(SEQ5=substr(Motif, 1, 5)) %>%
		dplyr::select(Category=Type, SEQ5, MU_5=ERV_rel_rate)

	rates3 <- ratelist[[2]] %>%
		mutate(SEQ3=substr(Motif, 1, 3)) %>%
		dplyr::select(Category=Type, SEQ3, MU_3=ERV_rel_rate)

	rates1 <- ratelist[[1]] %>%
		dplyr::select(Category=Type, MU_1=ERV_rel_rate)

	rates7A <- avrates %>%
		mutate(SEQ7=substr(Motif, 1, 7)) %>%
		dplyr::select(Category=Type, SEQ7, MU_7A=eur)

	chrp_c <- merge(chrp_c, rates7A, by=c("Category", "SEQ7"), all.x=T)
	chrp_c <- merge(chrp_c, rates7, by=c("Category", "SEQ7"), all.x=T)
	chrp_c <- merge(chrp_c, rates5, by=c("Category", "SEQ5"), all.x=T)
	chrp_c <- merge(chrp_c, rates3, by=c("Category", "SEQ3"))
	chrp_c <- merge(chrp_c, rates1, by=c("Category"))

	if(exists("maskgpdat")){
	rates_mask <- maskgpdat %>%
		mutate(SEQ7=substr(Motif, 1, 7)) %>%
		dplyr::select(Category=Type, SEQ7, MU_7M=ERV_rel_rate_mask)

		chrp_c <- merge(chrp_c, rates_mask, by=c("Category", "SEQ7"), all.x=T)
	}

	rates_7C <- r5m %>%
	mutate(SEQ7=substr(Motif, 1, 7)) %>%
	dplyr::select(Category=Type, SEQ7, MU_7C=MAC10_rel_rate)

	rates_7D <- r5m %>%
	mutate(SEQ7=substr(Motif, 1, 7)) %>%
	dplyr::select(Category=Type, SEQ7, MU_7D=ERV_down_rel_rate)

	# rates_anc <- ancgpdat %>%
	#   mutate(SEQ7=substr(Motif, 1, 7)) %>%
	#   dplyr::select(Category=Type, SEQ7, MU_7AN=ERV_rel_rate_anc)

	rates_anc <- test_anc %>%
	mutate(SEQ7=substr(Motif, 1, 7)) %>%
	dplyr::select(Category=Type, SEQ7, MU_7AN=estimate2)

	chrp_c <- merge(chrp_c, rates_7C, by=c("Category", "SEQ7"), all.x=T)
	chrp_c <- merge(chrp_c, rates_7D, by=c("Category", "SEQ7"), all.x=T)
	chrp_c <- merge(chrp_c, rates_anc, by=c("Category", "SEQ7"), all.x=T)

	chrp_c <- chrp_c %>%
		mutate(MU_7adj=-lambertW(-MU_7))

	chrp_c <- chrp_c %>%
	mutate(Category=ifelse(substr(SEQ,adj2+1,adj2+2)=="CG",
					paste0("cpg_",Category), Category)) %>%
	mutate(BIN=ceiling(POS/binw)) %>%
	# mutate(Category =
	#     plyr::mapvalues(Category, orderedcats1, orderedcats2)) %>%
	mutate(resid3=MU_3-MU_1,
		resid5=MU_5-MU_3,
		resid7=MU_7-MU_5,
		residL=MU-MU_7,
		resid5a=MU_5-MU_1,
		resid7a=MU_7-MU_1,
		residLa=MU-MU_1)

	return(chrp_c)
}

##############################################################################
# Pipeline:
# 1. sample non-mutated input sites
# 2. merge with DNMs
# 3. append columns for additional rate estimates
# 4. Simulate DNMs with same data
##############################################################################
validationPipe <- function(input_sites, input_dnms, nsites, ratelist, avrates){
	orderedcats <- c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA", "cpg_GC_AT", "cpg_GC_CG", "cpg_GC_TA")
	orderedcats1 <- c("AT_GC", "GC_AT", "cpg_GC_AT", "AT_CG", "GC_CG", "cpg_GC_CG", "AT_TA", "GC_TA", "cpg_GC_TA")
	orderedcats2 <- c("A>G", "C>T", "CpG>TpG", "A>C", "C>G", "CpG>GpG", "A>T", "C>A", "CpG>ApG")

	cat("Sampling non-mutated sites...\n")
	sampled_sites <- buildValidationData(input_sites, nsites)

	cat("Merging DNMs with sub-sampled background...\n")
	eval_sites <- bind_rows(list(sampled_sites, input_dnms)) %>%
		group_by(Category) %>%
		mutate(prop=cumsum(OBS)/sum(OBS)) %>%
		arrange(MU, prop) %>%
		mutate(SEQ7 = substr(SEQ, 1, 7),
			SEQ5 = substr(SEQ, 2, 6),
			SEQ3 = substr(SEQ, 3, 5))

	cat("Appending additional rate estimates...\n")
	eval_sites <- mergeRates(eval_sites, avrates, ratelist)

	cat("Generating simulated dataset...\n")
	simulated_dnms <- simMu(data=eval_sites, EST="MU", nobs=nrow(input_dnms), rseed=rseed)

	eval_sites_sim <- bind_rows(
		list(eval_sites[eval_sites$OBS==0,], simulated_dnms)) %>%
		group_by(Category) %>%
		mutate(prop=cumsum(OBS)/sum(OBS)) %>%
		arrange(MU, prop) %>%
		mutate(OBS=SIMOBS)

	##############################################################################
	# run overall models
	##############################################################################
	cat("Running combined models...\n")
	combined_models_sim <- runDNMLogit(eval_sites_sim, "COMBINED")
	combined_models <- runDNMLogit(eval_sites, "COMBINED")

	# Likelihood ratio test between nested models
	test13 <- lmtest::lrtest(combined_models[[1]], combined_models[[2]]) 
	test35 <- lmtest::lrtest(combined_models[[2]], combined_models[[3]])
	test57 <- lmtest::lrtest(combined_models[[3]], combined_models[[4]])
	test7L <- lmtest::lrtest(combined_models[[4]], combined_models[[5]])
	fulllist<-list(test13, test35, test57, test7L)

	pvals <- lapply(fulllist, function(x) x[[5]][2])

	combined_models_summary <- tibble(
			group="COMBINED",
			category="ALL",
			nsites=nsites,
			rsq=getNR2(combined_models),
			rsq_sim=getNR2(combined_models_sim),
			aic=getAIC(combined_models),
			aic_sim=getAIC(combined_models_sim)) %>%
		tibble::rownames_to_column("mod")

	fitfile <- paste0(parentdir, "/output/model_fit_", nsites, ".txt")
	write.table(combined_models_summary, fitfile, col.names=T, row.names=F, quote=F, sep="\t")

	##############################################################################
	# run type-specific models
	##############################################################################
	cat("Running type-specific models...\n")
	type_models_summary <- data.frame()
	lrtestdat <- data.frame()
	for(i in 1:length(orderedcats)){
		categ <- orderedcats2[i]

		type_dat_s <- eval_sites_sim %>%
			ungroup() %>%
			mutate(Category =
				plyr::mapvalues(Category, orderedcats1, orderedcats2)) %>%
			filter(Category==categ) %>%
			mutate(resid5=MU_5-MU_3, resid7=MU_7-MU_5, residL=MU-MU_7)
		type_models_sim <- runDNMLogit(type_dat_s, "TYPE")

		type_dat <- eval_sites %>%
			mutate(Category =
					plyr::mapvalues(Category, orderedcats1, orderedcats2)) %>%
			filter(Category==categ) %>%
			mutate(resid5=MU_5-MU_3, resid7=MU_7-MU_5, residL=MU-MU_7)
		type_models <- runDNMLogit(type_dat, "TYPE")

		test53 <- lmtest::lrtest(type_models[[1]], type_models[[2]])
		test75 <- lmtest::lrtest(type_models[[2]], type_models[[3]])
		testL7 <- lmtest::lrtest(type_models[[3]], type_models[[4]])
		fulllist <- list(test53, test75, testL7)
		pvals <- lapply(fulllist, function(x) x[[5]][2])

		mods <- c("5-mers", "7-mers", "Logit")
		lrtests <- tibble(category=categ, mod=mods, pvals)
		lrtestdat <- bind_rows(lrtestdat, lrtests)

		summ_tmp <- tibble(group="TYPE",
			category=categ,
			rsq=getNR2(type_models),
			rsq_sim=getNR2(type_models_sim),
			aic=getAIC(type_models),
				aic_sim=getAIC(type_models_sim)) %>%
			tibble::rownames_to_column("mod")

		type_models_summary <- bind_rows(type_models_summary, summ_tmp)
	}

	typefitfile <- paste0(parentdir, "/output/model_fit_by_type", nsites, ".txt")
	write.table(type_models_summary, typefitfile, col.names=T, row.names=F, quote=F, sep="\t")
	return(list(combined_models_summary, type_models_summary))
	}