####################
## Description:
##  - In this file we validation the semaglutide population
####################


## Load libraries
library(tidyverse)

## Set up directory path to save files (stagered to ensure folders are created)

dir.create("Samples")
dir.create("Samples/SGLT2-GLP1")

output_path <- "Samples/SGLT2-GLP1/Aurum"
dir.create(output_path)

## make directory for outputs
dir.create("Plots")

## male directory for outputs
dir.create(paste0(output_path, "/semaglutide"))


###############################################################################
###############################################################################
############################### Read Data In ##################################
###############################################################################
###############################################################################

## Load functions required

source("01.slade_aurum_functions.R")
source("02.slade_aurum_set_data.R")

###############################################################################
###############################################################################
########################## General variables ##################################
###############################################################################
###############################################################################

## Load dataset
semaglutide.dataset <- set_up_data_sglt2_glp1(dataset.type="semaglutide.dataset")

# load in variables used in the model
variables_mu <- readRDS(paste0(output_path, "/response_model_bcf/variables_mu.rds"))

variables_tau <- readRDS(paste0(output_path, "/response_model_bcf/variables_tau.rds"))

bcf_model <- readRDS(paste0(output_path, "/response_model_bcf/bcf_model.rds"))

# Create interim dataset with the variables needed
interim.dataset <- semaglutide.dataset %>%
  # select variables to make prediction
  select(all_of(c("patid", "pated", "ethnicity", "drugclass", "posthba1cfinal", unique(c(variables_mu, variables_tau))))) %>%
  drop_na()

# BCF has a quirk where you can only make predicted if both drugs are in the dataset
interim.dataset <- rbind(interim.dataset,
                         interim.dataset %>%
                           slice(1) %>%
                           mutate(drugclass = factor("SGLT2", levels = c("GLP1", "SGLT2"))))


# Predict treatment effect for these patients from our model
if (class(try(
  
  patient_effects <- readRDS(paste0(output_path, "/semaglutide/patient_effects.rds"))
  
  # predictions.interim <- readRDS(paste0(output_path, "/semaglutide/predictions.interim.rds"))
  
  , silent = TRUE)) == "try-error") {
  
  require(bcf)
  
  predictions.interim <- predict(object = bcf_model,
                                 x_predict_control = interim.dataset %>%
                                   select(
                                     all_of(variables_mu)
                                   ) %>%
                                   mutate_all(funs(as.numeric(.))) %>%
                                   as.matrix(),
                                 x_predict_moderate = interim.dataset %>%
                                   select(
                                     all_of(variables_tau)
                                   ) %>%
                                   mutate_all(funs(as.numeric(.))) %>%
                                   as.matrix(),
                                 pi_pred = rep(0.5, interim.dataset %>% nrow()),
                                 z_pred = interim.dataset %>%
                                   select(drugclass) %>%
                                   mutate(drugclass = ifelse(drugclass == "GLP1", 0, 1)) %>%
                                   unlist(),
                                 save_tree_directory = paste0(output_path, "/response_model_bcf/trees_no_prop"))
  
  
  saveRDS(predictions.interim, paste0(output_path, "/semaglutide/predictions.interim.rds"))
  
  # Combine and remove the extra iteration added
  patient_effects <- interim.dataset %>%
    select(patid, pated) %>%
    cbind(effects = colMeans(predictions.interim$tau)) %>%
    slice(-nrow(interim.dataset))
  
  
  saveRDS(patient_effects, paste0(output_path, "/semaglutide/patient_effects.rds"))
  
}



#:-------------------------------------------------------------------------------
#:-------------------------------------------------------------------------------
#:-------------------------------------------------------------------------------

### Validation of treatment effects - test to check if the model needs any adjustment (intercept or intercept + slope)

#closed testing function
closedtest <- function(cohort,dataset,observed,predicted,p.value){
  
  #Original model
  
  #Residuals
  resid <- observed-predicted
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #original log-likelihood
  n <- length(resid)
  logLik.original <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #Update intercept
  
  #Model with updated intercept
  m <- lm(observed-predicted~1,data=dataset)
  
  #Extract coefficient
  m1.intercept <- cbind(m$coefficients[1],confint(m)[1],confint(m)[2])
  
  #Residuals (actual - predicted)
  resid <- residuals(m)
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #log-likelihood
  n <- length(resid)
  logLik.intercept <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #Update slope & intercept
  
  #Model with updated slope &  intercept
  m <- lm(observed~predicted,data=dataset)
  
  #Extract coefficient
  m2.intercept <- cbind(m$coefficients[1],confint(m)[1,1],confint(m)[1,2])
  m2.slope <- cbind(m$coefficients[2],confint(m)[2,1],confint(m)[2,2])
  
  #Residuals (actual - predicted)
  resid <- residuals(m)
  
  #variance of residuals
  sigma2 <- var(resid)
  
  #log-likelihood
  n <- length(resid)
  logLik.recal <- -n/2 * log(2*pi*sigma2) - 1/(2*sigma2) * sum(resid^2)
  
  #test significance
  
  #1. Test recal in the large against the original model (no extra coeffs estimated)
  #2. If 2. is significant, test full recal against the recal in the large model using p df + 1 (1 extra coef estimated)
  #3. If 3. is significant, select 3. as final model, if not select 2. If neither, select 1 
  
  #ll diff
  dev_intercept <- -2*logLik.original + 2*logLik.intercept
  dev_recal <- -2*logLik.intercept + 2*logLik.recal
  
  #Diff in ll
  ncolx <- ncolx
  test1 <- (1-pchisq(dev_intercept, ncolx)) < p.value
  test2 <- (1-pchisq(dev_recal, ncolx+1)) < p.value
  
  #p.value
  p1 <- (1-pchisq(dev_intercept, ncolx))
  p2 <- (1-pchisq(dev_recal, ncolx+1))
  
  #Which model is chosen
  test_intercept <- 1 * (!test1)
  test_recal <- 2 * ((!test1)&(!test2))
  
  index_test <- (test_intercept + test_recal)
  
  res <- data.frame(cohort=c(cohort,cohort,cohort),
                    n=c(nrow(dataset),nrow(dataset),nrow(dataset)),
                    model=c("Original","Updated intercept","Recalibrated"),
                    loglikelihood=c(logLik.original,logLik.intercept,logLik.recal),
                    intercept=c(NA,
                                m1.intercept[1],
                                m1.intercept[2]),
                    intercept.w.ci=c(NA,
                                     paste0(round(m1.intercept[1],2), " (",paste0(round(m1.intercept[2],2),", ",paste0(round(m1.intercept[3],2)),")")),
                                     paste0(round(m2.intercept[1],2), " (",paste0(round(m2.intercept[2],2),", ",paste0(round(m2.intercept[3],2)),")"))),
                    
                    slope=c(NA,NA,m2.slope[1]),
                    slope.w.ci=c(NA,
                                 NA,
                                 paste0(round(m2.slope[1],2), " (",paste0(round(m2.slope[2],2),", ",paste0(round(m2.slope[3],2)),")"))),
                    p.value=c(NA,
                              round(p1,5),
                              round(p2,5)),
                    model.selected=c(ifelse(test1==FALSE & test2==FALSE,"Yes","No"),
                                     ifelse(test1==TRUE & test2==FALSE,"Yes","No"),
                                     ifelse(test2==TRUE,"Yes","No"))
  )
  return(res)
}


# global settings
p.value <- 0.05
ncolx <- 10
sample_frac <- 1

cohort <- "semaglutide"

observed <- interim.dataset$posthba1cfinal[-nrow(interim.dataset)]
predicted <- colMeans(predictions.interim$mu)[-nrow(interim.dataset)]

# dataset required for the calibration
dataset <- semaglutide.dataset %>%
  # select variables to make prediction
  select(all_of(c("patid", "pated", "ethnicity", "drugclass", "posthba1cfinal", unique(c(variables_mu, variables_tau))))) %>%
  drop_na()

# Run test
semaglutide.test <- closedtest(cohort,dataset,observed,predicted,p.value)

# Calculate the new adjusted predicted treatment effects
adjusted_effect <- patient_effects %>%
  mutate(effects = effects - semaglutide.test$intercept[2]) %>%
  mutate(best_drug = ifelse(effects > 0, "Favours GLP1", "Favours SGLT2"))

# Table of new optimal therapy
adjusted_effect_table <- adjusted_effect %>%
  select(best_drug) %>%
  table()/nrow(adjusted_effect)


# Plot new predicted treatment effects
adjusted_effect_hist <- hist_plot(adjusted_effect %>% rename("mean" = "effects"), 
                                  xmin = -10, xmax = 35, 
                                  title = paste0("Semaglutide cohort (n=", nrow(adjusted_effect), ") - adjusted intercept (SGLT2=", round(adjusted_effect_table[2]*100), "%, GLP1=", round(adjusted_effect_table[1]*100), "%)"))

# PDF with new predicted treatment effect histogram
pdf(width = 7, height = 5, "Plots/11.11.semaglutide_validation.pdf")
adjusted_effect_hist
dev.off()
