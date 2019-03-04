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
		NagelkerkeR2(x)))[seq(2, 2 * length(models), by = 2)]
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
mergeRates <- function(chrp_c, ratelist){
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
