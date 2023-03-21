####################
## Description:
##  - In this file we have the functions used for the analysis of Aurum
####################


new.vars <- function(data){
  # data = dataset that needs new vars
  
  data.new <- data %>%
    left_join(set_up_data_sglt2_glp1(dataset.type = "full.cohort"), by = c("patid", "pated")) %>%
    mutate(microvascular_complications = ifelse(prediabeticnephropathy == "Yes" | preneuropathy == "Yes" | preretinopathy == "Yes", "Yes", "No"),
           CV_problems = ifelse(preangina == "Yes" | preihd == "Yes" | premyocardialinfarction == "Yes" | prepad == "Yes" | prerevasc == "Yes" | prestroke == "Yes" | pretia == "Yes" | preaf == "Yes", "Yes", "No"),
           ASCVD = ifelse(premyocardialinfarction == "Yes" | prestroke == "Yes" | preihd == "Yes" | prepad == "Yes" | prerevasc == "Yes", "Yes", "No"),
           deprivation = factor(deprivation, labels = c("1", "1", "2", "2", "3", "3", "4", "4", "5", "5")),
           preckd = factor(preckd, labels = c("stage_1/2", "stage_1/2", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4", "stage_3a/stage_3b/stage_4")),
           drug_canagliflozin = ifelse(drugsubstances == "Canagliflozin" | drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Canagliflozin & Empagliflozin", "Yes", "No"),
           drug_dapagliflozin = ifelse(drugsubstances == "Canagliflozin & Dapagliflozin" | drugsubstances == "Dapagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin", "Yes", "No"),
           drug_empagliflozin = ifelse(drugsubstances == "Canagliflozin & Empagliflozin" | drugsubstances == "Dapagliflozin & Empagliflozin" | drugsubstances == "Empagliflozin", "Yes", "No"),
           drug_ertugliflozin = ifelse(drugsubstances == "Ertugliflozin", "Yes", "No"),
           drug_dulaglutide = ifelse(drugsubstances == "Dulaglutide" | drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Dulaglutide & Liraglutide", "Yes", "No"),
           drug_exenatide_short = ifelse(drugsubstances == "Dulaglutide & Exenatide" | drugsubstances == "Exenatide" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Liraglutide", "Yes", "No"),
           drug_exenatide_long = ifelse(drugsubstances == "Dulaglutide & Exenatide prolonged-release" | drugsubstances == "Exenatide & Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release" | drugsubstances == "Exenatide prolonged-release & Liraglutide", "Yes", "No"),
           drug_liraglutide = ifelse(drugsubstances == "Dulaglutide & Liraglutide" | drugsubstances == "Exenatide & Liraglutide" | drugsubstances == "Exenatide prolonged-release & Liraglutide" | drugsubstances == "Liraglutide" | drugsubstances == "Liraglutide & Lixisenatide", "Yes", "No"),
           drug_lixisenatide = ifelse(drugsubstances == "Liraglutide & Lixisenatide" | drugsubstances == "Lixisenatide", "Yes", "No"),
           prehba1c_na = ifelse(is.na(prehba1c), "Yes", "No"),
           prebmi_na = ifelse(is.na(prebmi), "Yes", "No"),
           preegfr_na = ifelse(is.na(preegfr), "Yes", "No"),
           prehdl_na = ifelse(is.na(prehdl), "Yes", "No"),
           prealt_na = ifelse(is.na(prealt), "Yes", "No"),
           prealbuminblood_na = ifelse(is.na(prealbuminblood), "Yes", "No"),
           prebilirubin_na = ifelse(is.na(prebilirubin), "Yes", "No"),
           pretotalcholesterol_na = ifelse(is.na(pretotalcholesterol), "Yes", "No"),
           premap_na = ifelse(is.na(premap), "Yes", "No"),
           posthba1cfinal_na = ifelse(is.na(posthba1cfinal), "Yes", "No"),
           hba1cmonth_na = ifelse(is.na(hba1cmonth), "Yes", "No"))
  
  return(data.new)
  
}





## Calculate assessments of prediction

rsq <- function (x, y) cor(x, y) ^ 2

calc_assessment <- function(data, posteriors, outcome_variable) {
  ##### Input variables
  # data - dataset used in the fitting 
  # posteriors - posteriors values for the dataset inputed
  # outcome_variable - variable with y values
  
  # Calculate R2
  r2 <- posteriors %>%
    apply(MARGIN = 1, function(x) rsq(data[,outcome_variable], x)) %>%
    quantile(probs = c(0.025, 0.5, 0.975))
  
  # Calculate RSS: residual sum of squares
  RSS <- posteriors %>%
    apply(MARGIN = 1, function(x) sum((data[,outcome_variable] - x)^2)) %>%
    quantile(probs = c(0.025, 0.5, 0.975))
  
  # Calculate RMSE: root mean square error
  RMSE <- posteriors %>%
    apply(MARGIN = 1, function(x) sqrt(sum((data[,outcome_variable] - x)^2)/nrow(data))) %>%
    quantile(probs = c(0.025, 0.5, 0.975))
  
  # return data.frame with all assessments
  assessment_values <- list(r2 = r2, RSS = RSS, RMSE = RMSE)
  
  return(assessment_values)
}


## Calculate residuals
calc_resid <- function(data, posteriors, outcome_variable) {
  ##### Input variables
  # data - dataset used in the fitting 
  # posteriors - posteriors values for the dataset inputed
  # outcome_variable - variable with outcome values
  
  # calculate standard deviation of residuals
  resid.SD <- apply(posteriors, MARGIN = 1, function(x) (data[,outcome_variable] - x)^2) %>%
    colSums() %>%
    as.data.frame() %>%
    set_names(c("SD")) %>%
    mutate(SD = sqrt(SD/(nrow(data)-2)))
  
  # calculate standardised residuals
  resid <- posteriors %>% t()
  for (i in 1:nrow(data)) {
    resid[i,] <- (data[i, outcome_variable] - resid[i,])/resid.SD[,1]
  }
  
  # return data.frame with residuals information for each data entry
  cred_pred <- cbind(lower_bd = apply(posteriors, MARGIN = 2, function(x) min(x)),
                     upper_bd = apply(posteriors, MARGIN = 2, function(x) max(x)),
                     mean = apply(posteriors, MARGIN = 2, function(x) mean(x)),
                     orig = data[,outcome_variable]) %>%
    as.data.frame() %>%
    mutate(resid = orig - mean,
           resid.low = orig - lower_bd,
           resid.high = orig - upper_bd) %>%
    cbind(std.resid = apply(resid, MARGIN = 1, function(x) mean(x)),
          std.resid.low = apply(resid, MARGIN = 1, function(x) min(x)),
          std.resid.high = apply(resid, MARGIN = 1, function(x) max(x)))
  
  return(cred_pred)
}


## Plot predicted vs observed and standardised residuals
resid_plot <- function(pred_dev, pred_val, title) {
  ##### Imput variables
  # pred_dev - predicted/observed values for development dataset
  # pred_val - predicted/observed values for validation dataset
  # title - plot title
  # 
  # # Plot of predicted vs observed for development dataset
  # plot_dev_pred <- pred_dev %>%
  #   ggplot() +
  #   theme_bw() +
  #   geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
  #   geom_point(aes(x = orig, y = mean)) +
  #   geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), linewidth=0.75) +
  #   xlim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
  #   ylim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
  #   xlab("Observed HbA1c (mmol/mol)") +
  #   ylab("Predicted HbA1c (mmol/mol)")
  # 
  # # Plot of predicted vs observed for validation dataset
  # plot_val_pred <- pred_val %>%
  #   ggplot() +
  #   theme_bw() +
  #   geom_errorbar(aes(ymin = lower_bd, ymax = upper_bd, x = orig), colour = "grey") +
  #   geom_point(aes(x = orig, y = mean)) +
  #   geom_abline(aes(intercept = 0, slope = 1), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), linewidth=0.75) +
  #   xlim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
  #   ylim(min(pred_dev$orig, pred_val$orig), max(pred_dev$orig, pred_val$orig)) +
  #   xlab("Observed HbA1c (mmol/mol)") +
  #   ylab("Predicted HbA1c (mmol/mol)")
  # 
  # Plot of standardised residuals for development dataset
  plot_dev_std <- pred_dev %>%
    ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = std.resid.low, ymax = std.resid.high, x = mean), colour = "grey") +
    geom_point(aes(x = mean, y = std.resid)) +
    geom_hline(aes(yintercept = 0), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), linewidth=0.75) +
    xlim(min(pred_dev$mean, pred_val$mean), max(pred_dev$mean, pred_val$mean)) +
    ylim(min(pred_dev$std.resid.low, pred_val$std.resid.low), max(pred_dev$std.resid.high, pred_val$std.resid.high)) +
    xlab("Average Predicted HbA1c (mmol/mol)") +
    ylab("Standardised Residuals")
  
  # Plot of standardised residuals for validation dataset
  plot_val_std <- pred_val %>%
    ggplot() +
    theme_bw() +
    geom_errorbar(aes(ymin = std.resid.low, ymax = std.resid.high, x = mean), colour = "grey") +
    geom_point(aes(x = mean, y = std.resid)) +
    geom_hline(aes(yintercept = 0), linetype ="dashed", color = viridis::viridis(1, begin = 0.6), linewidth=0.75) +
    xlim(min(pred_dev$mean, pred_val$mean), max(pred_dev$mean, pred_val$mean)) +
    ylim(min(pred_dev$std.resid.low, pred_val$std.resid.low), max(pred_dev$std.resid.high, pred_val$std.resid.high)) +
    xlab("Average Predicted HbA1c (mmol/mol)") +
    ylab("Standardised Residuals")
  
  plot_list <- list(plot_dev_std, plot_val_std)
  
  plot <- patchwork::wrap_plots(plot_list, ncol = 2) +
    patchwork::plot_annotation(
      tag_levels = "A", # labels A = development, B = validation
      title = title
    )
  
  return(plot)
  
}


hist_plot <- function(data, title, xmin, xmax, xtitle = "HbA1c difference (mmol/mol)", ytitle = "Number of people") {
  ### Input variables
  # data: dataset with column 'mean' corresponding to treatment effect
  # title: title for the plot
  # xmin: lower limit of x axis
  # xmax: upper limit of x axis
  # xtitle: title of x axis
  # ytitle: title of y axis
  
  # define data
  dat <- data %>% dplyr::select(mean) %>% mutate(above=ifelse(mean< 0, "Favours SGLT2i", "Favours GLP1-RA")) %>%
    mutate(above = factor(above, levels = c("Favours SGLT2i", "Favours GLP1-RA")))
  
  # plot
  plot <- ggplot(data = dat, aes(x = mean, fill = above)) +
    geom_histogram(position = "identity", alpha = 0.5, color = "black", breaks = seq(xmin, xmax, by = 1)) +
    geom_vline(aes(xintercept = 0), linetype = "dashed")+
    labs(title = title, x = xtitle, y = ytitle) +
    scale_fill_manual(values = c("#f1a340", "dodgerblue2"))+
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal")
  
  return(plot)
}



calc_ATE <- function(data, validation_type, variable, quantile_var, breakdown = NULL, prop_scores, caliper = 0.05, replace = FALSE, order = "random") {
  ##### Explanation of the function:
  # This function performs the calibration necessary for the model / clinical outcomes.
  # It can be done in three ways:
  #     - Propensity score matching:
  #           After splitting the population into subgroups, patients are matched
  #       based on their propensity scores. The matched patients are used to calculate
  #       the average treatment effect (ATE) and compared to the mean conditional
  #       average treatment effect (CATE).
  #     - Propensity score matching + adjustment:
  #           After splitting the population into subgroups, patients are matched
  #       based on their propensity scores. The matched patients are used to calculate
  #       the average treatment effect (ATE) whilst adjusting for all variables used
  #       in the model and compared to the mean conditional average treatment 
  #       effect (CATE).
  #     - Adjustment:
  #           After splitting the population into subgroups, The average treatment 
  #       effect (ATE) is calculated whilst adjusting for all variables used in the
  #       the model and it is compared to the mean conditional average treatment
  #       effect (CATE).
  
  ##### Input variables
  # data - data with variables + treatment effect quantiles 'quantile_var'
  # validation_type - string containing the type of validation you want to perform,
  #       that been: 
  #           - Propensity score matching: PSM 
  #           - Propensity score matching with variable adjustment: PSM + adjust
  #           - Variable adjustment: Adjust
  # variable - variable with outcome values
  # quantile_var - variable containing quantile/subgrouping indexes
  # breakdown - variables used to adjust  of matching
  # prop_scores - propensity scores for individuals or vector with variables from dataset
  # caliper - maximum distance between propensity scores of drug 1 vs drug 2
  # replace - logical variables, whether we replace matched individuals of small group
  # order - which side of propensity scores we start matching individuals, "largest", "smallest", "random"
  
  
  ##### Initial checks required for running the function:
  
  # If 'validation_type' is not supplied, error.
  if (missing(validation_type)) {stop("'validation_type' needs to be supplied")}
  # If 'validation_type' is not a character string, error.
  if (!is.character(validation_type)) {stop("'validation_type' must be a character string")}
  # If 'validation_type' is not one of the options in this list, error
  if (!(validation_type %in% c("PSM", "PSM + adjust", "Adjust"))) {
    stop("'validation_type' must be one of: 'PSM' / 'PSM + adjust' / 'Adjust'")
  }
  
  # If 'variable' is not supplied, error.
  if (missing(variable)) {stop("'variable' needs to be supplied")}
  # If 'variable' is not a character string, error.
  if (!is.character(variable)) {stop("'variable' must be a character string")}
  
  # If 'quantile_var' is not supplied, error.
  if (missing(quantile_var)) {stop("'quantile_var' needs to be supplied")}
  # If 'quantile_var' is not a character string, error.
  if (!is.character(quantile_var)) {stop("'quantile_var' must be a character string")}
  
  # load libraries
  require(rlang)
  
  ##### Start of the function
  
  # split predicted treatment effects into deciles
  predicted_treatment_effect <- data %>%
    plyr::ddply(quantile_var, dplyr::summarise,
                N = length(hba1c_diff),
                hba1c_diff.pred = mean(hba1c_diff))
  
  # maximum number of deciles being tested
  quantiles <- length(unique(data[,quantile_var]))
  
  # create lists with results
  mnumber = c(1:quantiles)
  models  <- as.list(1:quantiles)
  obs <- vector(); lci <- vector(); uci <- vector()
  
  if (validation_type == "PSM") {
    # This is the code for the calibration using PSM
    
    ##### Initial checks required for running the function:
    
    # If 'prop_scores' is not supplied, error.
    if (missing(prop_scores)) {stop("'prop_scores' needs to be supplied")}
    # If 'breakdown' is not supplied, error.
    if (is.null(breakdown)) {stop("'breakdown' needs to be supplied")}
    
    # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
    prop_score <- 1 - prop_scores
    
    # create lists with results
    n_drug1 <- vector(); n_drug2 <- vector(); matchit.ouputs <- list()
    
    # join propensity scores into dataset
    data.new <- data %>%
      cbind(prop_score)
    
    # iterate through deciles
    for (i in mnumber) {
      
      # model if propensity scores are provided
      matching_package_result <- MatchIt::matchit(
        formula = formula(paste0("drugclass ~ ", paste(breakdown, collapse = " + "))), # shouldn't be used since we are specifying 'distance' (propensity scores)
        data = data.new[which(data.new[,quantile_var] == i),], # select people in the quantile
        method = "nearest",
        distance = data.new[which(data.new[,quantile_var] == i),"prop_score"],
        replace = replace,
        m.order = order,
        caliper = caliper,
        mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
      )
      
      # collect matching summary
      matchit.ouputs[[i]] <- matching_package_result
      
      # calculate number of patients with either drug in the subgroup
      n_drug1 <- append(n_drug1, sum(!is.na(matching_package_result$match.matrix)))
      n_drug2 <- append(n_drug2, length(unique(matching_package_result$match.matrix[complete.cases(matching_package_result$match.matrix)])))
      
      # formula
      formula <- "posthba1cfinal ~ factor(drugclass)"
      
      # fit linear regression for decile in the matched dataset
      models[[i]] <- lm(as.formula(formula),data=data.new[data.new[,quantile_var] == i,], weights = matching_package_result$weights)
      
      # collect treatment effect from regression
      obs <- append(obs,models[[i]]$coefficients[2])
      
      # calculate confidence intervals
      confint_all <- confint(models[[i]], levels=0.95)
      
      # collect lower bound CI
      lci <- append(lci,confint_all[2,1])
      
      # collect upper bound CI
      uci <- append(uci,confint_all[2,2])
      
    }
    
    # join treatment effects for deciles in a data.frame
    effects <- data.frame(predicted_treatment_effect,cbind(n_drug1, n_drug2, obs, lci, uci))
    
    # returned list with fitted propensity model + decile treatment effects
    final_dataset <- list(effects = effects,
                          matching_outputs = matchit.ouputs)
    
    return(final_dataset)
    
  } else if (validation_type == "PSM + adjust") {
    # This is the code for the calibration using PSM + adjust
    
    ##### Initial checks required for running the function:
    
    # If 'prop_scores' is not supplied, error.
    if (missing(prop_scores)) {stop("'prop_scores' needs to be supplied")}
    # If 'breakdown' is not supplied, error.
    if (is.null(breakdown)) {stop("'breakdown' needs to be supplied")}
    
    # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
    prop_score <- 1 - prop_scores
    
    # create lists with results
    n_drug1 <- vector(); n_drug2 <- vector(); matchit.ouputs <- list()
    
    # join propensity scores into dataset
    data.new <- data %>%
      cbind(prop_score)
    
    # iterate through deciles
    for (i in mnumber) {
      
      # model if propensity scores are provided
      matching_package_result <- MatchIt::matchit(
        formula = formula(paste0("drugclass ~ ", paste(breakdown, collapse = " + "))), # shouldn't be used since we are specifying 'distance' (propensity scores)
        data = data.new[which(data.new[,quantile_var] == i),], # select people in the quantile
        method = "nearest",
        distance = data.new[which(data.new[,quantile_var] == i),"prop_score"],
        replace = replace,
        m.order = order,
        caliper = caliper,
        mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
      )
      
      # collect matching summary
      matchit.ouputs[[i]] <- matching_package_result
      
      # calculate number of patients with either drug in the subgroup
      n_drug1 <- append(n_drug1, sum(!is.na(matching_package_result$match.matrix)))
      n_drug2 <- append(n_drug2, length(unique(matching_package_result$match.matrix[complete.cases(matching_package_result$match.matrix)])))
      
      # variables used in adjustment
      breakdown_adjust <- breakdown
      # variables with more than one category represented
      checker <- which(sapply(data.new[data.new[,quantile_var] == i,breakdown_adjust], function(col) length(unique(col))) > 1)
      
      # formula
      formula <- paste0("posthba1cfinal ~ factor(drugclass) +", paste(breakdown_adjust[checker], collapse = "+"))
      
      # fit linear regression for decile in the matched dataset
      models[[i]] <- lm(as.formula(formula),data=data.new[data.new[,quantile_var] == i,], weights = matching_package_result$weights)
      
      # collect treatment effect from regression
      obs <- append(obs,models[[i]]$coefficients[2])
      
      # calculate confidence intervals
      confint_all <- confint(models[[i]], levels=0.95)
      
      # collect lower bound CI
      lci <- append(lci,confint_all[2,1])
      
      # collect upper bound CI
      uci <- append(uci,confint_all[2,2])
      
    }
    
    # join treatment effects for deciles in a data.frame
    effects <- data.frame(predicted_treatment_effect,cbind(n_drug1, n_drug2, obs, lci, uci))
    
    # returned list with fitted propensity model + decile treatment effects
    final_dataset <- list(effects = effects,
                          matching_outputs = matchit.ouputs)
    
    return(final_dataset)
    
  } else {
    # This is the code for the calibration using Adjust
    
    ##### Initial checks required for running the function:
    
    # If 'breakdown' is not supplied, error.
    if (is.null(breakdown)) {stop("'breakdown' needs to be supplied")}
    
    # iterate through deciles
    for (i in mnumber) {
      
      # do this differently if quantile_var is categorical
      if (is.factor(data[,quantile_var])) {
        # dataset being used in this quantile
        data.new <- data[data[,quantile_var] == levels(data[,quantile_var])[i],]
      } else {
        # dataset being used in this quantile
        data.new <- data[data[,quantile_var] == i,]
      }
      
      # variables used in adjustment
      breakdown_adjust <- breakdown
      # variables with only one variable represented
      checker <- which(sapply(data.new[,breakdown_adjust], function(col) length(unique(col))) > 1)
      
      # formula
      formula <- paste0("posthba1cfinal ~ factor(drugclass) +", paste(breakdown_adjust[checker], collapse = "+"))
      
      # fit linear regression for decile in the matched dataset
      models[[i]] <- lm(as.formula(formula),data=data.new)
      
      # collect treatment effect from regression
      obs <- append(obs,models[[i]]$coefficients[2])
      
      # calculate confidence intervals
      confint_all <- confint(models[[i]], levels=0.95)
      
      # collect lower bound CI
      lci <- append(lci,confint_all[2,1])
      
      # collect upper bound CI
      uci <- append(uci,confint_all[2,2])
      
    }
    
    # join treatment effects for deciles in a data.frame
    effects <- data.frame(predicted_treatment_effect,cbind(obs, lci, uci))
    
    # returned list with fitted propensity model + decile treatment effects
    final_dataset <- list(effects = effects)
    
    return(final_dataset)
    
  }
  
  
}















# calc_ATE_validation_adjust <- function(data, variable, quantile_var = "hba1c_diff.q", breakdown = NULL, adjust = TRUE) {
#   ##### Input variables
#   # data - Development dataset with variables + treatment effect quantiles (hba1c_diff.q)
#   # variable - variable with y values
#   # quantile_var - variable containing quantile indexes
#   # breakdown - variables used to compare quality of matching
#   # adjust - variables to adjust linear regression of effects
#   
#   require(rlang)
#   
#   # split predicted treatment effects into deciles
#   predicted_treatment_effect <- data %>%
#     plyr::ddply(quantile_var, dplyr::summarise,
#                 N = length(hba1c_diff),
#                 hba1c_diff.pred = mean(hba1c_diff))
#   
#   # maximum number of deciles being tested
#   quantiles <- length(unique(data[,quantile_var]))
#   
#   # create lists with results
#   mnumber = c(1:quantiles)
#   models  <- as.list(1:quantiles)
#   obs <- vector(); lci <- vector(); uci <- vector();
#   
#   
#   # iterate through deciles
#   for (i in mnumber) {
#     
#     # do this differently if quantile_var is categorical
#     if (is.factor(data[,quantile_var])) {
#       # dataset being used in this quantile
#       data.new <- data[data[,quantile_var] == levels(data[,quantile_var])[i],]
#     } else {
#       # dataset being used in this quantile
#       data.new <- data[data[,quantile_var] == i,]
#     }
#     
#     if (adjust == TRUE) {
#       
#       # variables used in adjustment
#       breakdown_adjust <- breakdown
#       # categorical variables in breakdown
#       # factors <- sapply(data.new[,breakdown_adjust], is.factor)
#       # variables with only one variable represented
#       checker <- which(sapply(data.new[,breakdown_adjust], function(col) length(unique(col))) > 1)
#       
#       formula <- paste0("posthba1cfinal ~ factor(drugclass) +", paste(breakdown_adjust[checker], collapse = "+"))
#       
#     } else {
#       formula <- "posthba1cfinal ~ factor(drugclass)"
#     }
#     
#     
#     
#     
#     
#     # fit linear regression for decile in the matched dataset
#     models[[i]] <- lm(as.formula(formula),data=data.new)
#     
#     # collect treatment effect from regression
#     obs <- append(obs,models[[i]]$coefficients[2])
#     
#     # calculate confidence intervals
#     confint_all <- confint(models[[i]], levels=0.95)
#     
#     # collect lower bound CI
#     lci <- append(lci,confint_all[2,1])
#     
#     # collect upper bound CI
#     uci <- append(uci,confint_all[2,2])
#     
#   }
#   
#   # join treatment effects for deciles in a data.frame
#   effects <- data.frame(predicted_treatment_effect,cbind(obs, lci, uci))
#   
#   # returned list with fitted propensity model + decile treatment effects
#   t <- list(effects = effects)
#   
#   return(t)
# }

# calc_ATE_validation_prop_matching <- function(data, variable, prop_scores, quantile_var="hba1c_diff.q", caliper = 0.05, replace = FALSE, order = "random", breakdown = NULL, adjust = FALSE) {
#   ##### Input variables
#   # data - Development dataset with variables + treatment effect quantiles (hba1c_diff.q)
#   # variable - variable with y values
#   # prop_scores - propensity scores for individuals or vector with variables from dataset
#   # quantile_var - variable containing quantile indexes
#   # caliper - maximum distance between propensity scores of drug 1 vs drug 2
#   # replace - logical variables, whether we replace matched individuals of small group
#   # order - which side we start matching individuals, "largest", "smallest", "random"
#   # breakdown - variables used to compare quality of matching
#   # adjust - variables to adjust linear regression of effects
#   
#   require(rlang)
#   
#   # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
#   prop_score <- 1 - prop_scores
#   
#   # split predicted treatment effects into deciles
#   predicted_treatment_effect <- data %>%
#     plyr::ddply(quantile_var, dplyr::summarise,
#                 N = length(hba1c_diff),
#                 hba1c_diff.pred = mean(hba1c_diff))
#   
#   # maximum number of deciles being tested
#   quantiles <- length(unique(data[,quantile_var]))
#   
#   # create lists with results
#   mnumber = c(1:quantiles)
#   models  <- as.list(1:quantiles)
#   obs <- vector(); lci <- vector(); uci <- vector(); n_drug1 <- vector(); n_drug2 <- vector()
#   
#   data.new <- data %>%
#     cbind(prop_score)
#   
#   # list of matchings
#   if (!is.null(breakdown)) {
#     matchit.ouputs <- list()
#   }
#   
#   # iterate through deciles
#   for (i in mnumber) {
#     
#     if (!is.null(breakdown)) {
#       # model if propensity scores are provided
#       matching_package_result <- MatchIt::matchit(
#         formula = formula(paste0("drugclass ~ ", paste(breakdown, collapse = " + "))), # shouldn't be used since we are specifying 'distance' (propensity scores)
#         data = data.new[which(data.new[,quantile_var] == i),], # select people in the quantile
#         method = "nearest",
#         distance = data.new[which(data.new[,quantile_var] == i),"prop_score"],
#         replace = replace,
#         m.order = order,
#         caliper = caliper,
#         mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
#       )
#     } else {
#       # model if propensity scores are provided
#       matching_package_result <- MatchIt::matchit(
#         formula = formula("drugclass ~ posthba1cfinal"), # shouldn't be used since we are specifying 'distance' (propensity scores)
#         data = data.new[which(data.new[,quantile_var] == i),], # select people in the quantile
#         method = "nearest",
#         distance = data.new[which(data.new[,quantile_var] == i),"prop_score"],
#         replace = replace,
#         m.order = order,
#         caliper = caliper,
#         mahvars = NULL, estimand = "ATT", exact = NULL, antiexact = NULL, discard = "none", reestimate = FALSE, s.weights = NULL, std.caliper = TRUE, ratio = 1, verbose = FALSE, include.obj = FALSE,
#       )
#     }
#     
#     if (!is.null(breakdown)) {
#       matchit.ouputs[[i]] <- matching_package_result
#     }
#     
#     n_drug1 <- append(n_drug1, sum(!is.na(matching_package_result$match.matrix)))
#     n_drug2 <- append(n_drug2, length(unique(matching_package_result$match.matrix[complete.cases(matching_package_result$match.matrix)])))
#     
#     if (!is.null(breakdown) & adjust == TRUE) {
#       
#       # variables used in adjustment
#       breakdown_adjust <- breakdown
#       # variables with only one variable represented
#       checker <- which(sapply(data.new[data.new[,quantile_var] == i,breakdown_adjust], function(col) length(unique(col))) > 1)
#       
#       
#       formula <- paste0("posthba1cfinal ~ factor(drugclass) +", paste(breakdown_adjust[checker], collapse = "+"))
#       
#     } else {
#       formula <- "posthba1cfinal ~ factor(drugclass)"
#     }
#     
#     # fit linear regression for decile in the matched dataset
#     models[[i]] <- lm(as.formula(formula),data=data.new[data.new[,quantile_var] == i,], weights = matching_package_result$weights)
#   
#     # collect treatment effect from regression
#     obs <- append(obs,models[[i]]$coefficients[2])
#     
#     # calculate confidence intervals
#     confint_all <- confint(models[[i]], levels=0.95)
#     
#     # collect lower bound CI
#     lci <- append(lci,confint_all[2,1])
#     
#     # collect upper bound CI
#     uci <- append(uci,confint_all[2,2])
#      
#   }
#     
#   
#   
#   # join treatment effects for deciles in a data.frame
#   effects <- data.frame(predicted_treatment_effect,cbind(n_drug1, n_drug2, obs, lci, uci))
#     
#   
#   if (!is.null(breakdown)) {
#     # returned list with fitted propensity model + decile treatment effects
#     t <- list(effects = effects,
#               matching_outputs = matchit.ouputs)
#   } else {
#     # returned list with fitted propensity model + decile treatment effects
#     t <- list(effects = effects)
#   }
#   
#   return(t)
# }

### inverse propensity score weighting 

calc_ATE_validation_inverse_prop_weighting <- function(data, variable, prop_scores, quantile_var="hba1c_diff.q") {
  ##### Input variables
  # data - Development dataset with variables + treatment effect quantiles (quantile_var)
  # variable - variable with y values
  # prop_scores - propensity scores for individuals
  # quantile_var - variable containing quantile indexes
  
  # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
  prop_score <- 1 - prop_scores
  
  # split predicted treatment effects into deciles
  predicted_treatment_effect <- data %>%
    plyr::ddply(quantile_var, dplyr::summarise,
                N = length(hba1c_diff),
                hba1c_diff.pred = mean(hba1c_diff))
  
  # maximum number of deciles being tested
  quantiles <- length(unique(data[,quantile_var]))
  
  # create lists with results
  mnumber = c(1:quantiles)
  models  <- as.list(1:quantiles)
  obs <- vector(); lci <- vector(); uci <- vector();
  
  # join dataset and propensity score
  data.new <- data %>%
    cbind(calc_prop = prop_score)
  
  # weights for SGLT2 Z = 1
  sglt2.data <- data.new %>%
    filter(drugclass == "SGLT2") %>%
    mutate(calc_prop = 1/(calc_prop))
  
  # weights for GLP1 Z = 0
  glp1.data <- data.new %>%
    filter(drugclass == "GLP1") %>%
    mutate(calc_prop = 1/(1-calc_prop))
  
  data.new <- rbind(sglt2.data, glp1.data)
  
  # formula
  formula <- paste0(variable, " ~ factor(drugclass)")
  
  # iterate through deciles
  for (i in mnumber) {
    # fit linear regression for decile
    models[[i]] <- lm(as.formula(formula),data=data.new,subset=data.new[,quantile_var]==i, weights = calc_prop)
    
    # collect treatment effect from regression
    obs <- append(obs,models[[i]]$coefficients[2])
    
    # calculate confidence intervals
    confint_all <- confint(models[[i]], levels=0.95)
    
    # collect lower bound CI
    lci <- append(lci,confint_all[2,1])
    
    # collect upper bound CI
    uci <- append(uci,confint_all[2,2])
    
  }
  
  # join treatment effects for deciles in a data.frame
  effects <- data.frame(predicted_treatment_effect,cbind(obs,lci,uci))
  
  # returned list with fitted propensity model + decile treatment effects  
  t <- list(effects = effects)
  
  return(t)
}

### inverse propensity score weighting stabilised

calc_ATE_validation_inverse_prop_weighting_stabilised <- function(data, variable, prop_scores, quantile_var="hba1c_diff.q") {
  ##### Input variables
  # data - Development dataset with variables + treatment effect quantiles (quantile_var)
  # variable - variable with y values
  # prop_scores - propensity scores for individuals
  # quantile_var - variable containing quantile indexes
  
  # keep propensity scores (1-score because bartMachine makes 1-GLP1 and 0-SGLT2, should be the way around)
  prop_score <- 1 - prop_scores
  
  # split predicted treatment effects into deciles
  predicted_treatment_effect <- data %>%
    plyr::ddply(quantile_var, dplyr::summarise,
                N = length(hba1c_diff),
                hba1c_diff.pred = mean(hba1c_diff))
  
  # maximum number of deciles being tested
  quantiles <- length(unique(data[,quantile_var]))
  
  # create lists with results
  mnumber = c(1:quantiles)
  models  <- as.list(1:quantiles)
  obs <- vector(); lci <- vector(); uci <- vector();
  
  # join dataset and propensity score
  data.new <- data %>%
    cbind(calc_prop = prop_score)
  
  # weights for SGLT2 Z = 1
  sglt2.data <- data.new %>%
    filter(drugclass == "SGLT2") %>%
    mutate(calc_prop = 1/(calc_prop))
  
  # stabilise propensity scores
  sglt2.data <- sglt2.data %>%
    mutate(calc_prop = calc_prop*(nrow(sglt2.data)/nrow(data.new)))
  
  
  # weights for GLP1 Z = 0
  glp1.data <- data.new %>%
    filter(drugclass == "GLP1") %>%
    mutate(calc_prop = 1/(1-calc_prop))
  
  # stabilise propensity scores
  glp1.data <- glp1.data %>%
    mutate(calc_prop = calc_prop*(nrow(glp1.data)/nrow(data.new)))
  
  
  data.new <- rbind(sglt2.data, glp1.data)
  
  # formula
  formula <- paste0(variable, " ~ factor(drugclass)")
  
  # iterate through deciles
  for (i in mnumber) {
    # fit linear regression for decile
    models[[i]] <- lm(as.formula(formula),data=data.new,subset=data.new[,quantile_var]==i, weights = calc_prop)
    
    # collect treatment effect from regression
    obs <- append(obs,models[[i]]$coefficients[2])
    
    # calculate confidence intervals
    confint_all <- confint(models[[i]], levels=0.95)
    
    # collect lower bound CI
    lci <- append(lci,confint_all[2,1])
    
    # collect upper bound CI
    uci <- append(uci,confint_all[2,2])
    
  }
  
  # join treatment effects for deciles in a data.frame
  effects <- data.frame(predicted_treatment_effect,cbind(obs,lci,uci))
  
  # returned list with fitted propensity model + decile treatment effects  
  t <- list(effects = effects)
  
  return(t)
}



#Function to output ATE by subgroup
ATE_plot <- function(data,pred,obs,obslowerci,obsupperci, ymin, ymax, colour_background = FALSE) {
  ###
  # data: dataset used in fitting,
  # pred: column with predicted values
  # obs: observed values
  # obslowerci: lower bound of CI for prediction
  # obsupperci: upper bound of CI for prediction
  # colour_background: colour of drug benefit on the background
  
  if (missing(ymin)) {
    ymin <- plyr::round_any(floor(min(c(unlist(data[obslowerci]), unlist(data[pred])))), 2, f = floor)
  }
  if (missing(ymax)) {
    ymax <- plyr::round_any(ceiling(max(c(unlist(data[obsupperci]), unlist(data[pred])))), 2, f = ceiling)
  }
  
  plot <- ggplot(data = data,aes_string(x = pred,y = obs))
  
  if (colour_background == TRUE) {
    
    plot <- plot +
      annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "dodgerblue2", alpha = 0.5)  + # top right
      annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= "#f1a340", alpha = 0.5) # bottom left
    
  }
  
  plot <- plot +
    geom_point(alpha = 1) + 
    theme_bw() +
    geom_errorbar(aes_string(ymin = obslowerci, ymax = obsupperci), colour = "black", width = 0.1) +
    ylab("Decile average treatment effect (mmol/mol)") + 
    xlab("Predicted conditional average treatment effect (mmol/mol)") +
    scale_x_continuous(limits = c(ymin, ymax), breaks = c(seq(ymin, ymax, by = 2))) +
    scale_y_continuous(limits = c(ymin, ymax), breaks = c(seq(ymin, ymax, by = 2))) +
    geom_abline(intercept = 0, slope = 1, color = "red", lwd = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") 
  
  return(plot)
  
}

# Function for grouping values into intervals
group_values <- function(data, variable, breaks) {
  ### Input variables
  # data: dataset used in splitting
  # variable: variable with values to be split
  # breaks: break points between values
  
  # stop in case 'variable' is not included in 'data'
  if (is.null(data[, variable])) {stop("'variable' not included in 'data'")}
  
  # include extra values so that extremes are included
  breaks.full <- c(breaks, floor(min(data[,variable], na.rm = TRUE)), ceiling(max(data[,variable], na.rm = TRUE)))
  
  new.data <- data %>%
    cbind(intervals = cut(data[, variable], breaks = breaks.full))
  
  return(new.data)
}

