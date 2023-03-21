####################
## Description:
##  - This file includes the framework for patient selection.
####################

# load libraries
require(tidyverse)


###############################################################################
###############################################################################
############################### Functions #####################################
###############################################################################
###############################################################################

set_up_data <- function(dataset.type, drugs = c("GLP1", "SGLT2")) {
  ##### Explanation of the function:
  # This function retrieve the original CPRD dataset, and applies all the exclusion
  # criteria required for obtaining one of the post-datasets required for analysis.
  # Throughout the function, there are snippets of code 'dataset.type == "diagnostics"'.
  # This provides a breakdown of the number of people excluded/collected at each stage.
  
  ##### Input variables
  # dataset.type: a character string mentioning the type of dataset required
  # drugs: a vector with the names of drugs required for the dataset
  
  ##### Initial checks required for running the function:
  
  # If 'dataset.type' is not supplied, error.
  if (missing(dataset.type)) {stop("'dataset.type' needs to be supplied")}
  # If 'dataset.type' is not a character string, error.
  if (!is.character(dataset.type)) {stop("'dataset.type' must be a character string")}
  # If 'dataset.type' is not one of the options in this list, error
  if (!(dataset.type %in% c("diagnostics", "synthetic", "full.cohort", "ps.model.train", "ps.model.test", "hba1c.train", "hba1c.test", "weight.dataset", "discontinuation.dataset", "egfr.dataset", "ckd.dataset" , "cvd.dataset", "hf.dataset", "no_co.dataset", "semaglutide.dataset", "micro_comp.dataset", "retinopathy.dataset", "insulin.dataset"))) {
    stop("'dataset.type' must be one of: diagnostics / synthetic / full.cohort / ps.model.train / ps.model.test / hba1c.train / hba1c.test / weight.dataset / discontinuation.dataset / egfr.dataset / ckd.dataset / cvd.dataset / hf.dataset / no_co.dataset / semaglutide.dataset / micro_comp.dataset / retinopathy.dataset / insulin.dataset")
  }
  # If 'drugs' is not supplied, error.
  if (missing(drugs)) {stop("'drugs' needs to be supplied")}
  # If 'drugs is not a character string, error.
  if (!is.character(drugs)) {stop("'drugs' must be a character string")}
  # If 'drugs' is not one of the options in this list, error.
  for (i in 1:length(drugs)) {
    if (!(drugs[i] %in% c("DPP4", "GLP1", "INS", "MFN", "SGLT2", "SU", "TZD"))) {
      stop("'drugs' must be one of: DPP4 / GLP1 / INS / MFN / SGLT2 / SU / TZD")
    }
  }
  
  ##### Start of the function:
  
  # load original dataset
  load("/slade/CPRD_data/mastermind_2022/20221205_t2d_1stinstance.Rda")
  
  cprd <- t2d_1stinstance
  
  
  ################################################
  ##### Select only input drugs
  ################################################
  
  cprd <- cprd %>% 
    filter(drugclass %in% drugs)
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Select only input drugs")
    print("################################################")
    print(nrow(cprd))
    print(table(cprd$drugclass))
    
  }
  
  ################################################
  ##### Drop patients initiating before 1/1/2013
  ################################################
  
  #######################
  # Explore adjusted HbA1c repsonse by calendar year
  
  cprd  <- cprd %>%
    mutate(yrdrugstart = format(dstartdate, format = "%Y")) %>%
    mutate(yrdrugstart = as.numeric(yrdrugstart))
  
  
  cprd <- cprd %>%
    mutate(dstartdate_cutoff = ifelse(dstartdate < "2013-01-01", 1 , NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop patients initiating before 1/1/2013")
    print("################################################")
    print(table(cprd$dstartdate_cutoff))
    print(table(cprd$dstartdate_cutoff, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(is.na(dstartdate_cutoff))
  
  
  ################################################
  ##### Drop if treated with insulin when starting new drug
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if treated with insulin when starting new drug")
    print("################################################")
    print(table(cprd$INS))
    print(table(cprd$INS, cprd$drugclass))
    
  }
  
  
  if (dataset.type == "insulin.dataset") {
    
    # select those treated with insulin
    cprd <- cprd %>% 
      filter(INS == 1)   
    
  } else {
    
    # select those not treated with insulin
    cprd <- cprd %>% 
      filter(INS == 0)      
    
  }
    
  ################################################
  ##### Drop patients with ESRD
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop patients with ESRD")
    print("################################################")
    print(table(cprd$preckdstage))
    print(table(cprd$preckdstage, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(preckdstage != "stage_5")
  
  
  ################################################
  ##### Drop if first-line treatment
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if first-line treatment")
    print("################################################")
    print(table(cprd$drugline_all))
    print(table(cprd$drugline_all, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(drugline_all != 1)
  
  
  ################################################
  ##### Drop if semaglutide
  ################################################
  
  cprd <- cprd %>%
    mutate(semaglutide_drug = ifelse(str_detect(drugsubstances, "Semaglutide"), 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if semaglutide")
    print("################################################")
    print(table(cprd$semaglutide_drug))
    
  }
  
  if (dataset.type == "semaglutide.dataset") {
    cprd <- cprd %>%
      filter(!is.na(semaglutide_drug))
  } else {
    cprd <- cprd %>%
      filter(is.na(semaglutide_drug))
  }
  
  ###############################################################################
  ###############################################################################
  ############################# Variable Prep ###################################
  ###############################################################################
  ###############################################################################
  
  ### Add variable that identifies an individual entry in the data
  
  cprd <- cprd %>%
    mutate(pated = paste(patid, drugclass, dstartdate, sep = ".")) %>%
    
    ################################################
  ##### Drug of interest
  ################################################
  
  mutate(drugclass = factor(drugclass, levels = drugs)) 
  
  ################################################
  ##### Outcome HbA1c # name: posthba1c12m (missing - 46383)
  #####   - posthba1c12m but if missing
  #####     - posthba1c6m
  ################################################
  
  cprd <- cprd %>%
    mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%
    mutate(posthba1cfinal = as.numeric(posthba1cfinal))
  
  
  ################################################
  ##### Sociodemographic variables
  ################################################
  
  cprd <- cprd %>%
    #####   - Age: agetx (new var)
    mutate(agetx = as.numeric(dstartdate_age)) %>%
    #####   - Sex: sex
    mutate(sex = factor(ifelse(gender == 1, "Male", "Female"))) %>%
    
    #####   - Duration of diabetes: t2dmduration
    mutate(t2dmduration = as.numeric(dstartdate_dm_dur_all)) %>%
    #####   - Ethnicity: ethnicity
    mutate(ethnicity = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "South Asian", "Black", "Other", "Mixed"))) %>%
    
    #####   - Deprivation: deprivation
    mutate(deprivation = factor(imd2015_10)) %>%
    
    #####   - Smoking Status: smoke
    mutate(smoke = factor(smoking_cat)) %>%
    
    #####   - Line Therapy: drugline: turn all > 4 to 5+
    mutate(drugline = ifelse(drugline_all > 4, 5, drugline_all)) %>%
    mutate(drugline = factor(drugline, levels = c(2, 3, 4, 5), labels = c("2", "3", "4", "5+"))) %>%
    
    #####   - Hospitalisations in previous year
    mutate(prehospitalisation = factor(hosp_admission_prev_year, levels = c(0, 1), labels = c("No", "Yes")))
  
  ################################################
  ##### Diabetes treatment
  ################################################
  
  cprd <- cprd %>%
    #####   - Drugs taken alongside treatment
    mutate(ncurrtx = DPP4 + SGLT2 + GLP1 + TZD + SU + MFN) %>%
    mutate(ncurrtx = ifelse(ncurrtx > 4, 5, ncurrtx)) %>%
    mutate(ncurrtx = factor(ncurrtx, levels = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5+"))) %>%
    
    #####   - Outcome month: hba1cmonth
    mutate(hba1cmonth_12 = difftime(posthba1c12mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth_6 = difftime(posthba1c6mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
    mutate(hba1cmonth = as.numeric(hba1cmonth))
  
  
  ################################################
  ##### Biomarkers
  ################################################
  
  #####   - hba1c: prehba1c (Nothing to do)
  #####   - BMI: prebmi (Nothing to do)
  #####   - eGFR: preegfr (Nothing to do)
  #####   - Albumin:Creatine ratio: preacr (Nothing to do)
  #####   - Serum albumin: prealbumin_blood (remove the _ from the name)
  
  cprd <- cprd %>%
    rename("prealbuminblood" = "prealbumin_blood",
           "prealbuminblooddate" = "prealbumin_blooddate",
           "prealbuminblooddrugdiff" = "prealbumin_blooddrugdiff")
  
  #####   - Alanine aminotransferase: prealt (Nothing to do)
  #####   - Aspartate aminotransferase: preast (Nothing to do)
  #####   - Bilirubin: prebilirubin (Nothing to do)
  #####   - Fasting glucose: prefastingglucose (Nothing to do)
  #####   - Fasting haematocrit: prehaematocrit (Nothing to do)
  #####   - Fasting haemoglobin: prehaemoglobin (Nothing to do)
  #####   - High-density lipoprotein (HDL): prehdl (Nothing to do)
  
  #####   - Mean arterial BP: premap
  
  cprd <- cprd %>%
    mutate(premap = predbp + ((presbp - predbp) / 3))
  
  #####   - Total cholesterol: pretotalcholesterol (Nothing to do)
  #####   - Triglycerides: pretriglyceride (Nothing to do)
  
  
  
  ################################################
  ##### Comorbidities
  ################################################
  
  cprd <- cprd %>%
    #####   - Angina: predrug_earliest_angina
    mutate(preangina = factor(predrug_angina, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Chronic Liver Disease: predrug_earliest_cld
    mutate(precld = factor(predrug_cld, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Diabetic Nephropathy: predrug_earliest_diabeticnephropathy
    mutate(prediabeticnephropathy = factor(predrug_diabeticnephropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Heart failure: predrug_earliest_heartfailure
    mutate(preheartfailure = factor(predrug_heartfailure, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Hypertension: predrug_earliest_hypertension
    mutate(prehypertension = factor(predrug_hypertension, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Ischaemic Heart Disease: predrug_earliest_ihd
    mutate(preihd = factor(predrug_ihd, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Myocardial Infarction: predrug_earliest_myocardialinfarction
    mutate(premyocardialinfarction = factor(predrug_myocardialinfarction, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Neuropathy: predrug_earliest_neuropathy
    mutate(preneuropathy = factor(predrug_neuropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Peripheral Arterial Disease: predrug_earliest_pad
    mutate(prepad = factor(predrug_pad, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Retinopathy: predrug_earliest_retinopathy
    mutate(preretinopathy = factor(predrug_retinopathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Cardiac Revascularisation: predrug_earliest_revasc
    mutate(prerevasc = factor(predrug_revasc, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Stroke: predrug_earliest_stroke
    mutate(prestroke = factor(predrug_stroke, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Transient Ischaemic Attack: predrug_earliest_tia
    mutate(pretia = factor(predrug_tia, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Atrial fibrillation: predrug_earliest_af
    mutate(preaf = factor(predrug_af, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - CKD stages
    mutate(preckd = factor(preckdstage)) %>%
    #####   - Pre-existing CVD
    mutate(predrug_cvd = ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1, 1, 0)) %>%
    mutate(predrug_cvd = factor(predrug_cvd, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Date of CV / HF death for later outcomes
    mutate(cv_death_date_any_cause=if_else(!is.na(death_date) & !is.na(cv_death_any_cause) & cv_death_any_cause==1, death_date, as.Date(NA)),
           cv_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(cv_death_primary_cause) & cv_death_primary_cause==1, death_date, as.Date(NA)),
           hf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(hf_death_any_cause) & hf_death_any_cause==1, death_date, as.Date(NA)),
           hf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(hf_death_primary_cause) & hf_death_primary_cause==1, death_date, as.Date(NA))) %>%
    #####   - Define 5 years post-drug start date for later censoring
    mutate(five_years_post_dstart=dstartdate+(365.25*5))
  
  
  
  ################################################
  ##### Add in later GLP1/SGLT2/TZD drug starts needed for censoring
  ################################################
  
  
  
  if (dataset.type == "ckd.dataset" | dataset.type == "cvd.dataset" | dataset.type == "hf.dataset" | dataset.type == "no_co.dataset" | dataset.type == "diagnostics" | dataset.type == "full.cohort" | dataset.type == "micro_comp.dataset" | dataset.type == "retinopathy.dataset") {
    
    # Add in later GLP1/SGLT2/TZD drug starts needed for censoring
    load("/slade/CPRD_data/mastermind_2022/20221205_t2d_all_drug_periods.Rda")
    
    later_sglt2 <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="SGLT2") %>%
                    select(patid, next_sglt2=dstartdate)), by="patid") %>%
      filter(next_sglt2>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE), .groups = "drop_last") %>%
      ungroup()
    
    
    later_glp1 <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="GLP1") %>%
                    select(patid, next_glp1=dstartdate)), by="patid") %>%
      filter(next_glp1>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_glp1_start=min(next_glp1, na.rm=TRUE), .groups = "drop_last") %>%
      ungroup()
    
    
    later_tzd <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="TZD") %>%
                    select(patid, next_tzd=dstartdate)), by="patid") %>%
      filter(next_tzd>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_tzd_start=min(next_tzd, na.rm=TRUE), .groups = "drop_last") %>%
      ungroup()
    
    
    # Load in all instances of ckd outcomes
    load("/slade/CPRD_data/mastermind_2022/20230125_ckd_outcomes_all.Rda")
    
    
    cprd <- cprd %>%
      left_join(later_sglt2, by=c("patid", "dstartdate")) %>%
      left_join(later_glp1, by=c("patid", "dstartdate")) %>%
      left_join(later_tzd, by=c("patid", "dstartdate")) %>%
      # only left_join first instances of outcomes
      left_join(ckd_outcomes, by=c("patid", "dstartdate", "drugclass"))
    
  }
  
  
  
  
  ###############################################################################
  ###############################################################################
  #################### Final dataset - all patients #############################
  ###############################################################################
  ###############################################################################
  #
  # Add all variables necessary for ALL analysis in the paper.
  #
  
  if (dataset.type == "ckd.dataset" | dataset.type == "cvd.dataset" | dataset.type == "hf.dataset" | dataset.type == "no_co.dataset" | dataset.type == "diagnostics" | dataset.type == "full.cohort" | dataset.type == "micro_comp.dataset" | dataset.type == "retinopathy.dataset") {
    
    final.dataset <- cprd %>%
      select(
        # information regarding patient
        patid, pated, multi_drug_start, timeprevcombo, drugsubstances,
        # response hba1c
        posthba1cfinal,
        # therapies of interest
        drugclass,
        # background
        MFN, DPP4, GLP1, SGLT2, SU, TZD,
        # Sociodemographic features
        agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
        # Diabetes treatment 
        drugline, ncurrtx, hba1cmonth, dstartdate, dstopdate, yrdrugstart,
        # Biomarkers
        prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
        prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
        # Comorbidities
        preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
        preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf, preckd,
        # Weight analysis
        preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate,
        # eGFR analysis
        postegfr12m, postegfr6m,
        # discontinuation
        stopdrug_6m_3mFU,
        # CKD
        preckd, predrug_cvd, postckdstage345date, egfr40_or_ckd5,
        # CVD
        predrug_cvd, postdrug_first_primary_incident_mi, postdrug_first_primary_incident_stroke, cv_death_date_primary_cause, 
        five_years_post_dstart, death_date, next_sglt2_start, next_tzd_start, gp_record_end, next_glp1_start,
        # HF
        postdrug_first_primary_hhf, hf_death_date_primary_cause,
        # No comorbidities
        postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, postdrug_first_heartfailure, 
        hf_death_date_any_cause, qrisk2_10yr_score,
        # Microvascular complications
        postdrug_first_diabeticnephropathy, postdrug_first_neuropathy, postdrug_first_retinopathy
      ) %>%
      as.data.frame()
    
  } else {
    
    final.dataset <- cprd %>%
      select(
        # information regarding patient
        patid, pated, multi_drug_start, timeprevcombo, drugsubstances,
        # response hba1c
        posthba1cfinal,
        # therapies of interest
        drugclass,
        # background
        MFN, DPP4, GLP1, SGLT2, SU, TZD,
        # Sociodemographic features
        agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
        # Diabetes treatment 
        drugline, ncurrtx, hba1cmonth, dstartdate, dstopdate, yrdrugstart,
        # Biomarkers
        prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
        prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
        # Comorbidities
        preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
        preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf, preckd,
        # Weight analysis
        preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate,
        # eGFR analysis
        postegfr12m, postegfr6m,
        # discontinuation
        stopdrug_6m_3mFU,
        # No comorbidities
        qrisk2_10yr_score
      ) %>%
      as.data.frame()
    
  }
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Final dataset - all patients")
    print("################################################")
    print(nrow(final.dataset))
    print(table(final.dataset$drugclass))
    
  }
  
  # if full cohort was requested
  if (dataset.type == "full.cohort" | dataset.type == "semaglutide.dataset" | dataset.type == "insulin.dataset") {
    return(final.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################ Synthetic dataset ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # Create synthetic dataset
  if (dataset.type == "synthetic") {
    # load package
    require(synthpop)
    
    set.seed(123)
    syn.dataset <- synthpop::syn(final.dataset %>%
                                   select(
                                     # response hba1c
                                     posthba1cfinal,
                                     # therapies of interest
                                     drugclass,
                                     # Sociodemographic features
                                     agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
                                     # Diabetes treatment 
                                     drugline, ncurrtx, hba1cmonth, yrdrugstart,
                                     # Biomarkers
                                     prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
                                     prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
                                     # Comorbidities
                                     preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
                                     preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
                                   ) %>%
                                   sample_n(520),
                                 print.flag = FALSE
    )
    
    return(syn.dataset$syn)
    
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################## Propensity score model ###############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  ps.model.dataset <- final.dataset
  
  #:----------------------------------------------------
  # Select variables needed
  
  ps.model.dataset <- ps.model.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
      # Diabetes treatment 
      drugline, ncurrtx, yrdrugstart,
      # Biomarkers
      prehba1c, prebmi, preegfr, 
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  # Training dataset
  set.seed(123)
  ps.model.dataset.train <- ps.model.dataset %>%
    group_by(drugclass) %>%
    sample_frac(.6) %>%
    ungroup() %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Propensity score training cohort")
    print("################################################")
    print(nrow(ps.model.dataset.train))
    print(table(ps.model.dataset.train$drugclass))
    
  }
  
  if (dataset.type == "ps.model.train") {
    return(ps.model.dataset.train)
  }
  
  
  # Testing dataset
  ps.model.dataset.test <- subset(ps.model.dataset, !(pated %in% ps.model.dataset.train$pated)) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Propensity score testing cohort")
    print("################################################")
    print(nrow(ps.model.dataset.test))
    print(table(ps.model.dataset.test$drugclass))
    
  }
  
  if (dataset.type == "ps.model.test") {
    return(ps.model.dataset.test)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################### HbA1c model ###################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model")
    print("################################################")
    
  }
  
  hba1c.model.dataset.train <- ps.model.dataset.train %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  hba1c.model.dataset.test <- ps.model.dataset.test %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$multi_drug_start))
    print(table(hba1c.model.dataset.train$multi_drug_start, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$multi_drug_start))
    print(table(hba1c.model.dataset.test$multi_drug_start, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(multi_drug_start == 0)
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$timeprevcombo_less61))
    print(table(hba1c.model.dataset.train$timeprevcombo_less61, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$timeprevcombo_less61))
    print(table(hba1c.model.dataset.test$timeprevcombo_less61, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(timeprevcombo_less61))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(timeprevcombo_less61))
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$hb_extreme_53))
    print(table(hba1c.model.dataset.train$hb_extreme_53, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$hb_extreme_53))
    print(table(hba1c.model.dataset.test$hb_extreme_53, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(hb_extreme_53))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print("Training Cohort")
    print(table(is.na(hba1c.model.dataset.train$prehba1c)))
    print(table(is.na(hba1c.model.dataset.train$prehba1c), hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(is.na(hba1c.model.dataset.test$prehba1c)))
    print(table(is.na(hba1c.model.dataset.test$prehba1c), hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(prehba1c))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if post HbA1c missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post HbA1c missing")
    print("################################################")
    print("Training Cohort")
    print(table(is.na(hba1c.model.dataset.train$posthba1cfinal)))
    print(table(is.na(hba1c.model.dataset.train$posthba1cfinal), hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(is.na(hba1c.model.dataset.test$posthba1cfinal)))
    print(table(is.na(hba1c.model.dataset.test$posthba1cfinal), hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(posthba1cfinal))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(posthba1cfinal))
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  # Training dataset
  final.hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model - Training cohort")
    print("################################################")
    print("Training Cohort")
    print(nrow(final.hba1c.model.dataset.train))
    print(table(final.hba1c.model.dataset.train$drugclass))
    
  }
  
  if (dataset.type == "hba1c.train") {
    return(final.hba1c.model.dataset.train)
  }
  
  # Testing dataset
  final.hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model - Testing cohort")
    print("################################################")
    print("Testing Cohort")
    print(nrow(final.hba1c.model.dataset.test))
    print(table(final.hba1c.model.dataset.test$drugclass))
    
  }
  
  if (dataset.type == "hba1c.test") {
    return(final.hba1c.model.dataset.test)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################ Weight population ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Weight model")
    print("################################################")
    
  }
  
  weight.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(weight.dataset$multi_drug_start))
    print(table(weight.dataset$multi_drug_start, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(weight.dataset$timeprevcombo_less61))
    print(table(weight.dataset$timeprevcombo_less61, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(weight.dataset$hb_extreme_53))
    print(table(weight.dataset$hb_extreme_53, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(weight.dataset$prehba1c)))
    print(table(is.na(weight.dataset$prehba1c), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if Weight is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if Weight is missing")
    print("################################################")
    print(table(is.na(weight.dataset$preweight)))
    print(table(is.na(weight.dataset$preweight), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(preweight))
  
  ################################################
  ##### Drop if post Weight is missing
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(postweight = ifelse(is.na(postweight12m), postweight6m, postweight12m))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post Weight is missing")
    print("################################################")
    print(table(is.na(weight.dataset$postweight)))
    print(table(is.na(weight.dataset$postweight), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(postweight))
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.weight.dataset <- weight.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # Weight analysis
      preweight, postweight
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Weight model - final")
    print("################################################")
    print(nrow(final.weight.dataset))
    print(table(final.weight.dataset$drugclass))
    
  }
  
  if (dataset.type == "weight.dataset") {
    return(final.weight.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################## Discontinuation population ###########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Discontinuation model")
    print("################################################")
    
  }
  
  discontinuation.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(discontinuation.dataset$multi_drug_start))
    print(table(discontinuation.dataset$multi_drug_start, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(discontinuation.dataset$timeprevcombo_less61))
    print(table(discontinuation.dataset$timeprevcombo_less61, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(discontinuation.dataset$hb_extreme_53))
    print(table(discontinuation.dataset$hb_extreme_53, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(discontinuation.dataset$prehba1c)))
    print(table(is.na(discontinuation.dataset$prehba1c), discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.discontinuation.dataset <- discontinuation.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # discontinuation
      stopdrug_6m_3mFU
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Discontinuation model - final")
    print("################################################")
    print(nrow(final.discontinuation.dataset))
    print(table(final.discontinuation.dataset$drugclass))
    
  }
  
  if (dataset.type == "discontinuation.dataset") {
    return(final.discontinuation.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################## eGFR population ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### eGFR model")
    print("################################################")
    
  }
  
  egfr.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(egfr.dataset$multi_drug_start))
    print(table(egfr.dataset$multi_drug_start, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(egfr.dataset$timeprevcombo_less61))
    print(table(egfr.dataset$timeprevcombo_less61, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(egfr.dataset$hb_extreme_53))
    print(table(egfr.dataset$hb_extreme_53, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$prehba1c)))
    print(table(is.na(egfr.dataset$prehba1c), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if baseline eGRF is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if baseline eGRF is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$preegfr)))
    print(table(is.na(egfr.dataset$preegfr), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(preegfr))
  
  
  ################################################
  ##### Drop if post eGRF is missing
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(postegfr = ifelse(is.na(postegfr12m), postegfr6m, postegfr12m))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post eGRF is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$postegfr)))
    print(table(is.na(egfr.dataset$postegfr), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(postegfr))
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.egfr.dataset <- egfr.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # eGFR analysis
      postegfr
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### eGFR model - final")
    print("################################################")
    print(nrow(final.egfr.dataset))
    print(table(final.egfr.dataset$drugclass))
    
  }
  
  if (dataset.type == "egfr.dataset") {
    return(final.egfr.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################### CKD outcome population ##############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CKD outcome model")
    print("################################################")
    
  }
  
  ckd.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(ckd.dataset$multi_drug_start))
    print(table(ckd.dataset$multi_drug_start, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(ckd.dataset$timeprevcombo_less61))
    print(table(ckd.dataset$timeprevcombo_less61, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(ckd.dataset$hb_extreme_53))
    print(table(ckd.dataset$hb_extreme_53, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(ckd.dataset$prehba1c)))
    print(table(is.na(ckd.dataset$prehba1c), ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if CKD
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(no_ckd = ifelse(!is.na(preckd) & (preckd=="stage_3a" | preckd=="stage_3b" | preckd=="stage_4"), 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if CKD")
    print("################################################")
    print(table(ckd.dataset$no_ckd))
    print(table(ckd.dataset$no_ckd, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(no_ckd))
  
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(ckd.dataset$TZD, ckd.dataset$drugclass))
    print("GLP1 treated")
    print(table(ckd.dataset$GLP1, ckd.dataset$drugclass))
    print("SGLT2 treated")
    print(table(ckd.dataset$SGLT2, ckd.dataset$drugclass))
    
  }
  
  
  ckd.dataset <- ckd.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  ckd.dataset <- ckd.dataset %>% 
    mutate(postdrug_stage345 = pmin(postckdstage345date, na.rm = TRUE),
           postdrug_egfr40_or_ckd5 = pmin(egfr40_or_ckd5, na.rm = TRUE)) %>%
    
    
  mutate(postdrug_stage345_censdate = if_else(drugclass=="GLP1",
                                              pmin(five_years_post_dstart,
                                                   death_date,
                                                   next_sglt2_start,
                                                   next_tzd_start,
                                                   gp_record_end,
                                                   postdrug_stage345, na.rm=TRUE),
                                              if_else(drugclass=="SGLT2",
                                                      pmin(five_years_post_dstart,
                                                           death_date,
                                                           next_glp1_start,
                                                           next_tzd_start,
                                                           gp_record_end,
                                                           postdrug_stage345, na.rm=TRUE),
                                                      as.Date(NA))),
         
         postdrug_egfr40_or_ckd5_censdate = if_else(drugclass=="GLP1",
                                                    pmin(five_years_post_dstart,
                                                         death_date,
                                                         next_sglt2_start,
                                                         next_tzd_start,
                                                         gp_record_end,
                                                         postdrug_egfr40_or_ckd5, na.rm=TRUE),
                                                    if_else(drugclass=="SGLT2",
                                                            pmin(five_years_post_dstart,
                                                                 death_date,
                                                                 next_glp1_start,
                                                                 next_tzd_start,
                                                                 gp_record_end,
                                                                 postdrug_egfr40_or_ckd5, na.rm=TRUE),
                                                            as.Date(NA))),
         
         
         postdrug_stage345_censvar=ifelse(!is.na(postdrug_stage345) & postdrug_stage345_censdate==postdrug_stage345, 1, 0),
         
         postdrug_egfr40_or_ckd5_censvar=ifelse(!is.na(postdrug_egfr40_or_ckd5) & postdrug_egfr40_or_ckd5_censdate==postdrug_egfr40_or_ckd5, 1, 0),
         
         postdrug_stage345_censtime_yrs=as.numeric(difftime(postdrug_stage345_censdate, dstartdate, unit="days"))/365.25,
         
         postdrug_egfr40_or_ckd5_censtime_yrs=as.numeric(difftime(postdrug_egfr40_or_ckd5_censdate, dstartdate, unit="days"))/365.25
         
         )
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.ckd.dataset <- ckd.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CKD model - final")
    print("################################################")
    print(nrow(final.ckd.dataset))
    print(table(final.ckd.dataset$drugclass))
    
  }
  
  if (dataset.type == "ckd.dataset") {
    return(final.ckd.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################### CVD outcome population ##############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD outcome model")
    print("################################################")
    
  }
  
  cvd.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(cvd.dataset$multi_drug_start))
    print(table(cvd.dataset$multi_drug_start, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  cvd.dataset <- cvd.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(cvd.dataset$timeprevcombo_less61))
    print(table(cvd.dataset$timeprevcombo_less61, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  cvd.dataset <- cvd.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(cvd.dataset$hb_extreme_53))
    print(table(cvd.dataset$hb_extreme_53, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(cvd.dataset$prehba1c)))
    print(table(is.na(cvd.dataset$prehba1c), cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if CVD
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if CVD")
    print("################################################")
    print(table(cvd.dataset$predrug_cvd))
    print(table(cvd.dataset$predrug_cvd, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(predrug_cvd == "No")
  
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(cvd.dataset$TZD, cvd.dataset$drugclass))
    print("GLP1 treated")
    print(table(cvd.dataset$GLP1, cvd.dataset$drugclass))
    print("SGLT2 treated")
    print(table(cvd.dataset$SGLT2, cvd.dataset$drugclass))
    
  }
  
  
  cvd.dataset <- cvd.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  
  ################################################
  ##### Outcome variables
  ################################################
  
  ## MACE: narrow definition (hospitalisation/death): narrow ('incident') MI/stroke HES codes (primary cause only) + CV death in ONS death (primary cause only)
  
  cvd.dataset <- cvd.dataset %>%
    
    mutate(postdrug_mace=pmin(postdrug_first_primary_incident_mi, postdrug_first_primary_incident_stroke, cv_death_date_primary_cause, na.rm=TRUE)) %>%
    
    mutate(postdrug_mace_censdate=if_else(drugclass=="GLP1",
                                          pmin(five_years_post_dstart,
                                               death_date,
                                               next_sglt2_start,
                                               next_tzd_start,
                                               gp_record_end,
                                               postdrug_mace, na.rm=TRUE),
                                          
                                          if_else(drugclass=="SGLT2",
                                                  pmin(five_years_post_dstart,
                                                       death_date,
                                                       next_glp1_start,
                                                       next_tzd_start,
                                                       gp_record_end,
                                                       postdrug_mace, na.rm=TRUE),
                                                  as.Date(NA))),
           
           postdrug_mace_censvar=ifelse(!is.na(postdrug_mace) & postdrug_mace_censdate==postdrug_mace, 1, 0),
           
           postdrug_mace_censtime_yrs=as.numeric(difftime(postdrug_mace_censdate, dstartdate, unit="days"))/365.25)
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.cvd.dataset <- cvd.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.cvd.dataset))
    print(table(final.cvd.dataset$drugclass))
    
  }
  
  if (dataset.type == "cvd.dataset") {
    return(final.cvd.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ##################### Heart Failure outcome population ########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Heart Failure outcome model")
    print("################################################")
    
  }
  
  hf.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(hf.dataset$multi_drug_start))
    print(table(hf.dataset$multi_drug_start, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  hf.dataset <- hf.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(hf.dataset$timeprevcombo_less61))
    print(table(hf.dataset$timeprevcombo_less61, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  hf.dataset <- hf.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(hf.dataset$hb_extreme_53))
    print(table(hf.dataset$hb_extreme_53, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(hf.dataset$prehba1c)))
    print(table(is.na(hf.dataset$prehba1c), hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  ################################################
  ##### Drop if heart failure
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if heart failure")
    print("################################################")
    print(table(hf.dataset$preheartfailure))
    print(table(hf.dataset$preheartfailure, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(preheartfailure == "No")
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(hf.dataset$TZD, hf.dataset$drugclass))
    print("GLP1 treated")
    print(table(hf.dataset$GLP1, hf.dataset$drugclass))
    print("SGLT2 treated")
    print(table(hf.dataset$SGLT2, hf.dataset$drugclass))
    
  }
  
  
  hf.dataset <- hf.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  hf.dataset <- hf.dataset %>%
    
    ## HF: narrow definition (hospitalisation/death): HF HES codes (primary cause only) + HF death in ONS death (primary cause only)
    
    
    mutate(postdrug_hf=pmin(postdrug_first_primary_hhf, hf_death_date_primary_cause, na.rm=TRUE)) %>%
    
    mutate(postdrug_hf_censdate=if_else(drugclass=="GLP1",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_tzd_start,
                                             gp_record_end,
                                             postdrug_hf, na.rm=TRUE),
                                        
                                        if_else(drugclass=="SGLT2",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_glp1_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_hf, na.rm=TRUE),
                                                as.Date(NA))),
           
           postdrug_hf_censvar=ifelse(!is.na(postdrug_hf) & postdrug_hf_censdate==postdrug_hf, 1, 0),
           
           postdrug_hf_censtime_yrs=as.numeric(difftime(postdrug_hf_censdate, dstartdate, unit="days"))/365.25)
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.hf.dataset <- hf.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.hf.dataset))
    print(table(final.hf.dataset$drugclass))
    
  }
  
  if (dataset.type == "hf.dataset") {
    return(final.hf.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ################# No comorbidities outcome population ########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### No comorbidities outcome model")
    print("################################################")
    
  }
  
  no_co.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(no_co.dataset$multi_drug_start))
    print(table(no_co.dataset$multi_drug_start, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(no_co.dataset$timeprevcombo_less61))
    print(table(no_co.dataset$timeprevcombo_less61, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(no_co.dataset$hb_extreme_53))
    print(table(no_co.dataset$hb_extreme_53, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(no_co.dataset$prehba1c)))
    print(table(is.na(no_co.dataset$prehba1c), no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  ################################################
  ##### Drop if comorbidity
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(comorbidities = ifelse(predrug_cvd == "No" & preheartfailure == "No" & (is.na(preckd) | (preckd!="stage_3a" & preckd!="stage_3b" & preckd!="stage_4")), NA_real_, 1))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if comorbidities")
    print("################################################")
    print(table(no_co.dataset$comorbidities, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(comorbidities))
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(no_co.dataset$TZD, no_co.dataset$drugclass))
    print("GLP1 treated")
    print(table(no_co.dataset$GLP1, no_co.dataset$drugclass))
    print("SGLT2 treated")
    print(table(no_co.dataset$SGLT2, no_co.dataset$drugclass))
    
  }
  
  
  no_co.dataset <- no_co.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome variables
  ################################################
  
  
  no_co.dataset <- no_co.dataset %>%
    
    
    mutate(postdrug_mace=pmin(postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, na.rm=TRUE),
           postdrug_hf=pmin(postdrug_first_heartfailure, hf_death_date_any_cause, na.rm=TRUE),
           postdrug_stage345 = pmin(postckdstage345date, na.rm = TRUE),
           postdrug_egfr40_or_ckd5 = pmin(egfr40_or_ckd5, na.rm = TRUE)) %>%
    
    
    # Mace: broad definition: MI/stroke GP codes + broad MI/stroke HES codes (any cause) + CV death in ONS death (any cause)
    
    
    mutate(postdrug_mace_censdate=if_else(drugclass=="GLP1",
                                          pmin(five_years_post_dstart,
                                               death_date,
                                               next_sglt2_start,
                                               next_tzd_start,
                                               gp_record_end,
                                               postdrug_mace, na.rm=TRUE),
                                          
                                          if_else(drugclass=="SGLT2",
                                                  pmin(five_years_post_dstart,
                                                       death_date,
                                                       next_glp1_start,
                                                       next_tzd_start,
                                                       gp_record_end,
                                                       postdrug_mace, na.rm=TRUE),
                                                  as.Date(NA))),
           
           postdrug_mace_censvar=ifelse(!is.na(postdrug_mace) & postdrug_mace_censdate==postdrug_mace, 1, 0),
           
           postdrug_mace_censtime_yrs=as.numeric(difftime(postdrug_mace_censdate, dstartdate, unit="days"))/365.25,
           
           
           ## HF: broad definition: HF GP codes + HF HES codes (any cause) + HF death in ONS death (any cause)
           
           
           postdrug_hf_censdate=if_else(drugclass=="GLP1",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_tzd_start,
                                             gp_record_end,
                                             postdrug_hf, na.rm=TRUE),
                                        
                                        if_else(drugclass=="SGLT2",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_glp1_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_hf, na.rm=TRUE),
                                                as.Date(NA))),
           
           postdrug_hf_censvar=ifelse(!is.na(postdrug_hf) & postdrug_hf_censdate==postdrug_hf, 1, 0),
           
           postdrug_hf_censtime_yrs=as.numeric(difftime(postdrug_hf_censdate, dstartdate, unit="days"))/365.25,
           
           
           
           ## CKD: reaching stage 3/4/5 OR decrease of 40% in egfr and stage 5
           
           
           
           postdrug_stage345_censdate = if_else(drugclass=="GLP1",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_sglt2_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_stage345, na.rm=TRUE),
                                                if_else(drugclass=="SGLT2",
                                                        pmin(five_years_post_dstart,
                                                             death_date,
                                                             next_glp1_start,
                                                             next_tzd_start,
                                                             gp_record_end,
                                                             postdrug_stage345, na.rm=TRUE),
                                                        as.Date(NA))),
           
           postdrug_egfr40_or_ckd5_censdate = if_else(drugclass=="GLP1",
                                                      pmin(five_years_post_dstart,
                                                           death_date,
                                                           next_sglt2_start,
                                                           next_tzd_start,
                                                           gp_record_end,
                                                           postdrug_egfr40_or_ckd5, na.rm=TRUE),
                                                      if_else(drugclass=="SGLT2",
                                                              pmin(five_years_post_dstart,
                                                                   death_date,
                                                                   next_glp1_start,
                                                                   next_tzd_start,
                                                                   gp_record_end,
                                                                   postdrug_egfr40_or_ckd5, na.rm=TRUE),
                                                              as.Date(NA))),
           
           
           postdrug_stage345_censvar=ifelse(!is.na(postdrug_stage345) & postdrug_stage345_censdate==postdrug_stage345, 1, 0),
           
           postdrug_egfr40_or_ckd5_censvar=ifelse(!is.na(postdrug_egfr40_or_ckd5) & postdrug_egfr40_or_ckd5_censdate==postdrug_egfr40_or_ckd5, 1, 0),
           
           postdrug_stage345_censtime_yrs=as.numeric(difftime(postdrug_stage345_censdate, dstartdate, unit="days"))/365.25,
           
           postdrug_egfr40_or_ckd5_censtime_yrs=as.numeric(difftime(postdrug_egfr40_or_ckd5_censdate, dstartdate, unit="days"))/365.25
           
    )
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.no_co.dataset <- no_co.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.no_co.dataset))
    print(table(final.no_co.dataset$drugclass))
    
  }
  
  if (dataset.type == "no_co.dataset") {
    return(final.no_co.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ########## No microvascular complications outcome population ##################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### No microvascular complications outcome model")
    print("################################################")
    
  }
  
  micro_comp.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(micro_comp.dataset$multi_drug_start))
    print(table(micro_comp.dataset$multi_drug_start, micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  micro_comp.dataset <- micro_comp.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(micro_comp.dataset$timeprevcombo_less61))
    print(table(micro_comp.dataset$timeprevcombo_less61, micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  micro_comp.dataset <- micro_comp.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(micro_comp.dataset$hb_extreme_53))
    print(table(micro_comp.dataset$hb_extreme_53, micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(micro_comp.dataset$prehba1c)))
    print(table(is.na(micro_comp.dataset$prehba1c), micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  ################################################
  ##### Drop if comorbidity
  ################################################
  
  micro_comp.dataset <- micro_comp.dataset %>%
    mutate(comorbidities = ifelse(prediabeticnephropathy == "No" & preneuropathy == "No" & preretinopathy == "No", NA_real_, 1))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if comorbidities")
    print("################################################")
    print(table(micro_comp.dataset$comorbidities, micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(is.na(comorbidities))
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(micro_comp.dataset$TZD, micro_comp.dataset$drugclass))
    print("GLP1 treated")
    print(table(micro_comp.dataset$GLP1, micro_comp.dataset$drugclass))
    print("SGLT2 treated")
    print(table(micro_comp.dataset$SGLT2, micro_comp.dataset$drugclass))
    
  }
  
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  micro_comp.dataset <- micro_comp.dataset %>%
    
    ## Microvascular complications: diabetic nephropathy, retinopathy, neuropahty
    
    
    mutate(postdrug_micro_comp=pmin(postdrug_first_diabeticnephropathy, postdrug_first_neuropathy, postdrug_first_retinopathy, na.rm=TRUE)) %>%
    
    mutate(postdrug_micro_comp_censdate=if_else(drugclass=="GLP1",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_tzd_start,
                                             gp_record_end,
                                             postdrug_micro_comp, na.rm=TRUE),
                                        
                                        if_else(drugclass=="SGLT2",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_glp1_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_micro_comp, na.rm=TRUE),
                                                as.Date(NA))),
           
           postdrug_micro_comp_censvar=ifelse(!is.na(postdrug_micro_comp) & postdrug_micro_comp_censdate==postdrug_micro_comp, 1, 0),
           
           postdrug_micro_comp_censtime_yrs=as.numeric(difftime(postdrug_micro_comp_censdate, dstartdate, unit="days"))/365.25)
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.micro_comp.dataset <- micro_comp.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Microvascular complications model - final")
    print("################################################")
    print(nrow(final.micro_comp.dataset))
    print(table(final.micro_comp.dataset$drugclass))
    
  }
  
  if (dataset.type == "micro_comp.dataset") {
    return(final.micro_comp.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  #################### Retinopathy outcome population ###########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Retinopathy outcome model")
    print("################################################")
    
  }
  
  retinopathy.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(retinopathy.dataset$multi_drug_start))
    print(table(retinopathy.dataset$multi_drug_start, retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  retinopathy.dataset <- retinopathy.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(retinopathy.dataset$timeprevcombo_less61))
    print(table(retinopathy.dataset$timeprevcombo_less61, retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  retinopathy.dataset <- retinopathy.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(retinopathy.dataset$hb_extreme_53))
    print(table(retinopathy.dataset$hb_extreme_53, retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(retinopathy.dataset$prehba1c)))
    print(table(is.na(retinopathy.dataset$prehba1c), retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if retinopathy
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if retinopathy")
    print("################################################")
    print(table(retinopathy.dataset$preretinopathy))
    print(table(retinopathy.dataset$preretinopathy, retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(preretinopathy == "No")
  
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(retinopathy.dataset$TZD, retinopathy.dataset$drugclass))
    print("GLP1 treated")
    print(table(retinopathy.dataset$GLP1, retinopathy.dataset$drugclass))
    print("SGLT2 treated")
    print(table(retinopathy.dataset$SGLT2, retinopathy.dataset$drugclass))
    
  }
  
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  retinopathy.dataset <- retinopathy.dataset %>% 
    mutate(postdrug_retinopathy = pmin(postdrug_first_retinopathy, na.rm = TRUE)) %>%
    
    
    mutate(postdrug_retinopathy_censdate = if_else(drugclass=="GLP1",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_sglt2_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_retinopathy, na.rm=TRUE),
                                                if_else(drugclass=="SGLT2",
                                                        pmin(five_years_post_dstart,
                                                             death_date,
                                                             next_glp1_start,
                                                             next_tzd_start,
                                                             gp_record_end,
                                                             postdrug_retinopathy, na.rm=TRUE),
                                                        as.Date(NA))),
           
           postdrug_retinopathy_censvar=ifelse(!is.na(postdrug_retinopathy) & postdrug_retinopathy_censdate==postdrug_retinopathy, 1, 0),
           
           postdrug_retinopathy_censtime_yrs=as.numeric(difftime(postdrug_retinopathy_censdate, dstartdate, unit="days"))/365.25
           
    )
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.retinopathy.dataset <- retinopathy.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Retinopathy outcome model - final")
    print("################################################")
    print(nrow(final.retinopathy.dataset))
    print(table(final.retinopathy.dataset$drugclass))
    
  }
  
  if (dataset.type == "retinopathy.dataset") {
    return(final.retinopathy.dataset)
  }
  
  
  
  
  
}


set_up_data_sglt2_glp1 <- function(dataset.type) {
  ##### Explanation of the function:
  # This function retrieves the original CPRD dataset, and applies all the exclusion
  # criteria required for obtaining one of the post-datasets required for analysis.
  # Throughout the function, there are snippets of code 'dataset.type == "diagnositcs"'.
  # These provide a breakdown of the number of people excluded/collected at each stage.
  
  
  ##### Input variables
  # dataset.type: a character string mentioning the type of dataset required
  
  ##### Initial checks required for running the function:
  
  # If 'dataset.type' is not supplied, error.
  if (missing(dataset.type)) {stop("'dataset.type' needs to be supplied")}
  # If 'dataset.type' is not a character string, error.
  if (!is.character(dataset.type)) {stop("'dataset.type' must be a character string")}
  # If 'dataset.type' is not one of the options in this list, error.
  if (!(dataset.type %in% c("diagnostics", "synthetic", "full.cohort", "ps.model.train", "ps.model.test", "hba1c.train", "hba1c.test", "weight.dataset", "discontinuation.dataset", "egfr.dataset", "ckd.dataset" , "cvd.dataset", "hf.dataset", "no_co.dataset", "semaglutide.dataset", "micro_comp.dataset", "retinopathy.dataset", "insulin.dataset"))) {
    stop("'dataset.type' must be one of: diagnostics / synthetic / full.cohort / ps.model.train / ps.model.test / hba1c.train / hba1c.test / weight.dataset / discontinuation.dataset / egfr.dataset / ckd.dataset / cvd.dataset / hf.dataset / no_co.dataset / semaglutide.dataset / micro_comp.dataset / retinopathy.dataset / insulin.dataset")
  }
  
  ##### Start of the function:
  
  # load original dataset
  load("/slade/CPRD_data/mastermind_2022/20221205_t2d_1stinstance.Rda")
  
  cprd <- t2d_1stinstance
  
  
  ################################################
  ##### Select only SGLT-2 and GLP-1
  ################################################
  
  #SGLT2 vs GLP1
  cprd <- cprd %>% 
    filter(drugclass == "GLP1" | drugclass == "SGLT2")
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Select only SGLT-2 and GLP-1")
    print("################################################")
    print(nrow(cprd))
    print(table(cprd$drugclass))
    
  }
  
  
  ################################################
  ##### Drop patients initiating before 1/1/2013
  ################################################
  
  #######################
  # Explore adjusted HbA1c repsonse by calendar year
  
  cprd  <- cprd %>%
    mutate(yrdrugstart = format(dstartdate, format = "%Y")) %>%
    mutate(yrdrugstart = as.numeric(yrdrugstart))
  
  #######################
  
  
  cprd <- cprd %>%
    mutate(dstartdate_cutoff = ifelse(dstartdate < "2013-01-01", 1 , NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop patients initiating before 1/1/2013")
    print("################################################")
    print(table(cprd$dstartdate_cutoff))
    print(table(cprd$dstartdate_cutoff, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(is.na(dstartdate_cutoff))

  
  ################################################
  ##### Drop if treated with insulin when starting new drug
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if treated with insulin when starting new drug")
    print("################################################")
    print(table(cprd$INS))
    print(table(cprd$INS, cprd$drugclass))
    
  }
  
  if (dataset.type == "insulin.dataset") {
    
    # select those treated with insulin
    cprd <- cprd %>%
      filter(INS == 1)
    
  } else {
    
    # select those not treated with insulin
    cprd <- cprd %>% 
      filter(INS == 0)  
    
  }
  
  ################################################
  ##### Drop patients with ESRD
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop patients with ESRD")
    print("################################################")
    print(table(cprd$preckdstage))
    print(table(cprd$preckdstage, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(preckdstage != "stage_5")
  
  
  ################################################
  ##### Drop if first-line treatment
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if first-line treatment")
    print("################################################")
    print(table(cprd$drugline_all))
    print(table(cprd$drugline_all, cprd$drugclass))
    
  }
  
  cprd <- cprd %>%
    filter(drugline_all != 1)
  
  
  ################################################
  ##### Drop if semaglutide
  ################################################
  
  cprd <- cprd %>%
    mutate(semaglutide_drug = ifelse(str_detect(drugsubstances, "Semaglutide"), 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if semaglutide")
    print("################################################")
    print(table(cprd$semaglutide_drug))
    
  }
  
  if (dataset.type == "semaglutide.dataset") {
    cprd <- cprd %>%
      filter(!is.na(semaglutide_drug))
  } else {
    cprd <- cprd %>%
      filter(is.na(semaglutide_drug))
  }
  
  ###############################################################################
  ###############################################################################
  ############################# Variable Prep ###################################
  ###############################################################################
  ###############################################################################
  
  ### Add variable that identifies an individual entry in the data
  
  cprd <- cprd %>%
    mutate(pated = paste(patid, drugclass, dstartdate, sep = ".")) %>%
  
  ################################################
  ##### Drug of interest
  ################################################
  
    #####   - GLP1 and SGLT2: drugclass
    mutate(drugclass = factor(drugclass)) 
  # 1 - GLP1; 2 - SGLT2
  
  ################################################
  ##### Outcome HbA1c # name: posthba1c12m (missing - 46383)
  #####   - posthba1c12m but if missing
  #####     - posthba1c6m
  ################################################
  
  cprd <- cprd %>%
    mutate(posthba1cfinal = ifelse(is.na(posthba1c12m), posthba1c6m, posthba1c12m)) %>%
    mutate(posthba1cfinal = as.numeric(posthba1cfinal))
    
  
  ################################################
  ##### Sociodemographic variables
  ################################################
  
  cprd <- cprd %>%
    #####   - Age: agetx (new var)
    mutate(agetx = as.numeric(dstartdate_age)) %>%
    #####   - Sex: sex
    mutate(sex = factor(ifelse(gender == 1, "Male", "Female"))) %>%
    
    #####   - Duration of diabetes: t2dmduration
    mutate(t2dmduration = as.numeric(dstartdate_dm_dur_all)) %>%
    #####   - Ethnicity: ethnicity
    mutate(ethnicity = factor(ethnicity_5cat, levels = c(0, 1, 2, 3, 4), labels = c("White", "South Asian", "Black", "Other", "Mixed"))) %>%
    
    #####   - Deprivation: deprivation
    mutate(deprivation = factor(imd2015_10)) %>%
    
    #####   - Smoking Status: smoke
    mutate(smoke = factor(smoking_cat)) %>%
    
    #####   - Line Therapy: drugline: turn all > 4 to 5+
    mutate(drugline = ifelse(drugline_all > 4, 5, drugline_all)) %>%
    mutate(drugline = factor(drugline, levels = c(2, 3, 4, 5), labels = c("2", "3", "4", "5+"))) %>%
    
    #####   - Hospitalisations in previous year
    mutate(prehospitalisation = factor(hosp_admission_prev_year, levels = c(0, 1), labels = c("No", "Yes")))
  
  ################################################
  ##### Diabetes treatment
  ################################################
  
  cprd <- cprd %>%
    #####   - Drugs taken alongside treatment
    mutate(ncurrtx = DPP4 + SGLT2 + GLP1 + TZD + SU + MFN) %>%
    mutate(ncurrtx = ifelse(ncurrtx > 4, 5, ncurrtx)) %>%
    mutate(ncurrtx = factor(ncurrtx, levels = c(1, 2, 3, 4, 5), labels = c("1", "2", "3", "4", "5+"))) %>%
    
  #####   - Outcome month: hba1cmonth
    mutate(hba1cmonth_12 = difftime(posthba1c12mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth_6 = difftime(posthba1c6mdate, dstartdate, units = "days") / 30) %>%
    mutate(hba1cmonth = ifelse(is.na(hba1cmonth_12), hba1cmonth_6, hba1cmonth_12)) %>%
    mutate(hba1cmonth = as.numeric(hba1cmonth))
  
  
  ################################################
  ##### Biomarkers
  ################################################
  
    #####   - hba1c: prehba1c (Nothing to do)
    #####   - BMI: prebmi (Nothing to do)
    #####   - eGFR: preegfr (Nothing to do)
    #####   - Albumin:Creatine ratio: preacr (Nothing to do)
    #####   - Serum albumin: prealbumin_blood (remove the _ from the name)
  
  cprd <- cprd %>%
    rename("prealbuminblood" = "prealbumin_blood",
           "prealbuminblooddate" = "prealbumin_blooddate",
           "prealbuminblooddrugdiff" = "prealbumin_blooddrugdiff")
  
    #####   - Alanine aminotransferase: prealt (Nothing to do)
    #####   - Aspartate aminotransferase: preast (Nothing to do)
    #####   - Bilirubin: prebilirubin (Nothing to do)
    #####   - Fasting glucose: prefastingglucose (Nothing to do)
    #####   - Fasting haematocrit: prehaematocrit (Nothing to do)
    #####   - Fasting haemoglobin: prehaemoglobin (Nothing to do)
    #####   - High-density lipoprotein (HDL): prehdl (Nothing to do)
  
    #####   - Mean arterial BP: premap
  
  cprd <- cprd %>%
    mutate(premap = predbp + ((presbp - predbp) / 3))
  
    #####   - Total cholesterol: pretotalcholesterol (Nothing to do)
    #####   - Triglycerides: pretriglyceride (Nothing to do)
  
  
  
  ################################################
  ##### Comorbidities
  ################################################
  
  cprd <- cprd %>%
    #####   - Angina: predrug_earliest_angina
    mutate(preangina = factor(predrug_angina, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Chronic Liver Disease: predrug_earliest_cld
    mutate(precld = factor(predrug_cld, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Diabetic Nephropathy: predrug_earliest_diabeticnephropathy
    mutate(prediabeticnephropathy = factor(predrug_diabeticnephropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Heart failure: predrug_earliest_heartfailure
    mutate(preheartfailure = factor(predrug_heartfailure, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Hypertension: predrug_earliest_hypertension
    mutate(prehypertension = factor(predrug_hypertension, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Ischaemic Heart Disease: predrug_earliest_ihd
    mutate(preihd = factor(predrug_ihd, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Myocardial Infarction: predrug_earliest_myocardialinfarction
    mutate(premyocardialinfarction = factor(predrug_myocardialinfarction, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Neuropathy: predrug_earliest_neuropathy
    mutate(preneuropathy = factor(predrug_neuropathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Peripheral Arterial Disease: predrug_earliest_pad
    mutate(prepad = factor(predrug_pad, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Retinopathy: predrug_earliest_retinopathy
    mutate(preretinopathy = factor(predrug_retinopathy, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Cardiac Revascularisation: predrug_earliest_revasc
    mutate(prerevasc = factor(predrug_revasc, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Stroke: predrug_earliest_stroke
    mutate(prestroke = factor(predrug_stroke, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Transient Ischaemic Attack: predrug_earliest_tia
    mutate(pretia = factor(predrug_tia, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Atrial fibrillation: predrug_earliest_af
    mutate(preaf = factor(predrug_af, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - CKD stages
    mutate(preckd = factor(preckdstage)) %>%
    #####   - Pre-existing CVD
    mutate(predrug_cvd = ifelse(predrug_angina==1 | predrug_ihd==1 | predrug_myocardialinfarction==1 | predrug_pad==1 | predrug_revasc==1 | predrug_stroke==1, 1, 0)) %>%
    mutate(predrug_cvd = factor(predrug_cvd, levels = c(0, 1), labels = c("No", "Yes"))) %>%
    #####   - Date of CV / HF death for later outcomes
    mutate(cv_death_date_any_cause=if_else(!is.na(death_date) & !is.na(cv_death_any_cause) & cv_death_any_cause==1, death_date, as.Date(NA)),
           cv_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(cv_death_primary_cause) & cv_death_primary_cause==1, death_date, as.Date(NA)),
           hf_death_date_any_cause=if_else(!is.na(death_date) & !is.na(hf_death_any_cause) & hf_death_any_cause==1, death_date, as.Date(NA)),
           hf_death_date_primary_cause=if_else(!is.na(death_date) & !is.na(hf_death_primary_cause) & hf_death_primary_cause==1, death_date, as.Date(NA))) %>%
    #####   - Define 5 years post-drug start date for later censoring
    mutate(five_years_post_dstart=dstartdate+(365.25*5))
  
  
  
  ################################################
  ##### Add in later GLP1/SGLT2/TZD drug starts needed for censoring
  ################################################
  
  
  
  if (dataset.type == "ckd.dataset" | dataset.type == "cvd.dataset" | dataset.type == "hf.dataset" | dataset.type == "no_co.dataset" | dataset.type == "diagnostics" | dataset.type == "full.cohort" | dataset.type == "micro_comp.dataset" | dataset.type == "retinopathy.dataset") {
    
    # Add in later GLP1/SGLT2/TZD drug starts needed for censoring
    load("/slade/CPRD_data/mastermind_2022/20221205_t2d_all_drug_periods.Rda")
    
    later_sglt2 <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="SGLT2") %>%
                    select(patid, next_sglt2=dstartdate)), by="patid") %>%
      filter(next_sglt2>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_sglt2_start=min(next_sglt2, na.rm=TRUE), .groups = "drop_last") %>%
      ungroup()
    
    
    later_glp1 <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="GLP1") %>%
                    select(patid, next_glp1=dstartdate)), by="patid") %>%
      filter(next_glp1>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_glp1_start=min(next_glp1, na.rm=TRUE), .groups = "drop_last") %>%
      ungroup()
    
    
    later_tzd <- cprd %>%
      select(patid, dstartdate) %>%
      inner_join((t2d_all_drug_periods %>%
                    filter(drugclass=="TZD") %>%
                    select(patid, next_tzd=dstartdate)), by="patid") %>%
      filter(next_tzd>dstartdate) %>%
      group_by(patid, dstartdate) %>%
      summarise(next_tzd_start=min(next_tzd, na.rm=TRUE), .groups = "drop_last") %>%
      ungroup()
    
    
    # Load in all instances of ckd outcomes
    load("/slade/CPRD_data/mastermind_2022/20230125_ckd_outcomes_all.Rda")
    
    cprd <- cprd %>%
      left_join(later_sglt2, by=c("patid", "dstartdate")) %>%
      left_join(later_glp1, by=c("patid", "dstartdate")) %>%
      left_join(later_tzd, by=c("patid", "dstartdate")) %>%
      # only left_join first instances of outcomes
      left_join(ckd_outcomes, by=c("patid", "dstartdate", "drugclass"))
    
    
  }
  
  
  ###############################################################################
  ###############################################################################
  #################### Final dataset - all patients #############################
  ###############################################################################
  ###############################################################################
  #
  # Add all variables necessary for ALL analysis in the paper.
  #
  
  if (dataset.type == "ckd.dataset" | dataset.type == "cvd.dataset" | dataset.type == "hf.dataset" | dataset.type == "no_co.dataset" | dataset.type == "diagnostics" | dataset.type == "full.cohort" | dataset.type == "micro_comp.dataset" | dataset.type == "retinopathy.dataset") {
    
    final.dataset <- cprd %>%
      select(
        # information regarding patient
        patid, pated, multi_drug_start, timeprevcombo, drugsubstances,
        # response hba1c
        posthba1cfinal,
        # therapies of interest
        drugclass,
        # background
        MFN, DPP4, GLP1, SGLT2, SU, TZD,
        # Sociodemographic features
        agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
        # Diabetes treatment 
        drugline, ncurrtx, hba1cmonth, dstartdate, dstopdate, yrdrugstart,
        # Biomarkers
        prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
        prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
        # Comorbidities
        preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
        preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf, preckd,
        # Weight analysis
        preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate,
        # eGFR analysis
        postegfr12m, postegfr6m,
        # discontinuation
        stopdrug_6m_3mFU,
        # CKD
        preckd, predrug_cvd, postckdstage345date, egfr40_or_ckd5,
        # CVD
        predrug_cvd, postdrug_first_primary_incident_mi, postdrug_first_primary_incident_stroke, cv_death_date_primary_cause, 
        five_years_post_dstart, death_date, next_sglt2_start, next_tzd_start, gp_record_end, next_glp1_start,
        # HF
        postdrug_first_primary_hhf, hf_death_date_primary_cause,
        # No comorbidities
        postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, postdrug_first_heartfailure, 
        hf_death_date_any_cause, qrisk2_10yr_score,
        # Microvascular complications
        postdrug_first_diabeticnephropathy, postdrug_first_neuropathy, postdrug_first_retinopathy
      ) %>%
      as.data.frame()
    
  } else {
    
    final.dataset <- cprd %>%
      select(
        # information regarding patient
        patid, pated, multi_drug_start, timeprevcombo, drugsubstances,
        # response hba1c
        posthba1cfinal,
        # therapies of interest
        drugclass,
        # background
        MFN, DPP4, GLP1, SGLT2, SU, TZD,
        # Sociodemographic features
        agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
        # Diabetes treatment 
        drugline, ncurrtx, hba1cmonth, dstartdate, dstopdate, yrdrugstart,
        # Biomarkers
        prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
        prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
        # Comorbidities
        preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
        preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf, preckd,
        # Weight analysis
        preweight, postweight12m, postweight6m, postweight12mdate, postweight6mdate,
        # eGFR analysis
        postegfr12m, postegfr6m,
        # discontinuation
        stopdrug_6m_3mFU,
        # No comorbidities
        qrisk2_10yr_score
      ) %>%
      as.data.frame()
    
  }
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Final dataset - all patients")
    print("################################################")
    print(nrow(final.dataset))
    print(table(final.dataset$drugclass))
    
  }
  
  # if full cohort was requested
  if (dataset.type == "full.cohort" | dataset.type == "semaglutide.dataset" | dataset.type == "insulin.dataset") {
    return(final.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################ Synthetic dataset ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # Create synthetic dataset
  if (dataset.type == "synthetic") {
    # load package
    require(synthpop)
    
    set.seed(123)
    syn.dataset <- synthpop::syn(final.dataset %>%
                                   select(
                                     # response hba1c
                                     posthba1cfinal,
                                     # therapies of interest
                                     drugclass,
                                     # Sociodemographic features
                                     agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
                                     # Diabetes treatment 
                                     drugline, ncurrtx, hba1cmonth, yrdrugstart,
                                     # Biomarkers
                                     prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin, prefastingglucose,
                                     prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
                                     # Comorbidities
                                     preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
                                     preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
                                   ) %>%
                                   sample_n(520),
                                 print.flag = FALSE
    )
    
    return(syn.dataset$syn)
    
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################## Propensity score model ###############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  ps.model.dataset <- final.dataset
  
  #:----------------------------------------------------
  # Select variables needed
  
  ps.model.dataset <- ps.model.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, ethnicity, deprivation, smoke, prehospitalisation,
      # Diabetes treatment 
      drugline, ncurrtx, yrdrugstart,
      # Biomarkers
      prehba1c, prebmi, preegfr, 
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    )
  
  # Training dataset
  set.seed(123)
  ps.model.dataset.train <- ps.model.dataset %>%
    group_by(drugclass) %>%
    sample_frac(.6) %>%
    ungroup() %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Propensity score training cohort")
    print("################################################")
    print(nrow(ps.model.dataset.train))
    print(table(ps.model.dataset.train$drugclass))
    
  }
  
  if (dataset.type == "ps.model.train") {
    return(ps.model.dataset.train)
  }
  
  
  # Testing dataset
  ps.model.dataset.test <- subset(ps.model.dataset, !(pated %in% ps.model.dataset.train$pated)) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Propensity score testing cohort")
    print("################################################")
    print(nrow(ps.model.dataset.test))
    print(table(ps.model.dataset.test$drugclass))
    
  }
  
  if (dataset.type == "ps.model.test") {
    return(ps.model.dataset.test)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################### HbA1c model ###################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model")
    print("################################################")
    
  }
  
  hba1c.model.dataset.train <- ps.model.dataset.train %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  hba1c.model.dataset.test <- ps.model.dataset.test %>%
    select(patid, pated) %>%
    left_join(final.dataset, by = c("patid", "pated"))
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$multi_drug_start))
    print(table(hba1c.model.dataset.train$multi_drug_start, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$multi_drug_start))
    print(table(hba1c.model.dataset.test$multi_drug_start, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(multi_drug_start == 0)
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$timeprevcombo_less61))
    print(table(hba1c.model.dataset.train$timeprevcombo_less61, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$timeprevcombo_less61))
    print(table(hba1c.model.dataset.test$timeprevcombo_less61, hba1c.model.dataset.test$drugclass))
    
  }

  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(timeprevcombo_less61))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(timeprevcombo_less61))
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print("Training Cohort")
    print(table(hba1c.model.dataset.train$hb_extreme_53))
    print(table(hba1c.model.dataset.train$hb_extreme_53, hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(hba1c.model.dataset.test$hb_extreme_53))
    print(table(hba1c.model.dataset.test$hb_extreme_53, hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(is.na(hb_extreme_53))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print("Training Cohort")
    print(table(is.na(hba1c.model.dataset.train$prehba1c)))
    print(table(is.na(hba1c.model.dataset.train$prehba1c), hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(is.na(hba1c.model.dataset.test$prehba1c)))
    print(table(is.na(hba1c.model.dataset.test$prehba1c), hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(prehba1c))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if post HbA1c missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post HbA1c missing")
    print("################################################")
    print("Training Cohort")
    print(table(is.na(hba1c.model.dataset.train$posthba1cfinal)))
    print(table(is.na(hba1c.model.dataset.train$posthba1cfinal), hba1c.model.dataset.train$drugclass))
    print("Testing Cohort")
    print(table(is.na(hba1c.model.dataset.test$posthba1cfinal)))
    print(table(is.na(hba1c.model.dataset.test$posthba1cfinal), hba1c.model.dataset.test$drugclass))
    
  }
  
  hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    filter(!is.na(posthba1cfinal))
  
  hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    filter(!is.na(posthba1cfinal))
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  # Training dataset
  final.hba1c.model.dataset.train <- hba1c.model.dataset.train %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model - Training cohort")
    print("################################################")
    print("Training Cohort")
    print(nrow(final.hba1c.model.dataset.train))
    print(table(final.hba1c.model.dataset.train$drugclass))
    
  }
  
  if (dataset.type == "hba1c.train") {
    return(final.hba1c.model.dataset.train)
  }
  
  # Testing dataset
  final.hba1c.model.dataset.test <- hba1c.model.dataset.test %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration,
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### HbA1c model - Testing cohort")
    print("################################################")
    print("Testing Cohort")
    print(nrow(final.hba1c.model.dataset.test))
    print(table(final.hba1c.model.dataset.test$drugclass))
    
  }
  
  if (dataset.type == "hba1c.test") {
    return(final.hba1c.model.dataset.test)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################ Weight population ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Weight model")
    print("################################################")
    
  }
  
  weight.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(weight.dataset$multi_drug_start))
    print(table(weight.dataset$multi_drug_start, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(weight.dataset$timeprevcombo_less61))
    print(table(weight.dataset$timeprevcombo_less61, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(weight.dataset$hb_extreme_53))
    print(table(weight.dataset$hb_extreme_53, weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(weight.dataset$prehba1c)))
    print(table(is.na(weight.dataset$prehba1c), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if Weight is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if Weight is missing")
    print("################################################")
    print(table(is.na(weight.dataset$preweight)))
    print(table(is.na(weight.dataset$preweight), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(preweight))
  
  ################################################
  ##### Drop if post Weight is missing
  ################################################
  
  weight.dataset <- weight.dataset %>%
    mutate(postweight = ifelse(is.na(postweight12m), postweight6m, postweight12m))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post Weight is missing")
    print("################################################")
    print(table(is.na(weight.dataset$postweight)))
    print(table(is.na(weight.dataset$postweight), weight.dataset$drugclass))
    
  }
  
  weight.dataset <- weight.dataset %>%
    filter(!is.na(postweight))
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.weight.dataset <- weight.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # Weight analysis
      preweight, postweight
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Weight model - final")
    print("################################################")
    print(nrow(final.weight.dataset))
    print(table(final.weight.dataset$drugclass))
    
  }
  
  if (dataset.type == "weight.dataset") {
    return(final.weight.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################## Discontinuation population ###########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Discontinuation model")
    print("################################################")
    
  }
  
  discontinuation.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(discontinuation.dataset$multi_drug_start))
    print(table(discontinuation.dataset$multi_drug_start, discontinuation.dataset$drugclass))
    
  }

  discontinuation.dataset <- discontinuation.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(discontinuation.dataset$timeprevcombo_less61))
    print(table(discontinuation.dataset$timeprevcombo_less61, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  discontinuation.dataset <- discontinuation.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(discontinuation.dataset$hb_extreme_53))
    print(table(discontinuation.dataset$hb_extreme_53, discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(discontinuation.dataset$prehba1c)))
    print(table(is.na(discontinuation.dataset$prehba1c), discontinuation.dataset$drugclass))
    
  }
  
  discontinuation.dataset <- discontinuation.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.discontinuation.dataset <- discontinuation.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # discontinuation
      stopdrug_6m_3mFU
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Discontinuation model - final")
    print("################################################")
    print(nrow(final.discontinuation.dataset))
    print(table(final.discontinuation.dataset$drugclass))
    
  }
  
  if (dataset.type == "discontinuation.dataset") {
    return(final.discontinuation.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ############################## eGFR population ################################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### eGFR model")
    print("################################################")
    
  }
  
  egfr.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(egfr.dataset$multi_drug_start))
    print(table(egfr.dataset$multi_drug_start, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(egfr.dataset$timeprevcombo_less61))
    print(table(egfr.dataset$timeprevcombo_less61, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(egfr.dataset$hb_extreme_53))
    print(table(egfr.dataset$hb_extreme_53, egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$prehba1c)))
    print(table(is.na(egfr.dataset$prehba1c), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if baseline eGRF is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if baseline eGRF is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$preegfr)))
    print(table(is.na(egfr.dataset$preegfr), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(preegfr))
  
  
  ################################################
  ##### Drop if post eGRF is missing
  ################################################
  
  egfr.dataset <- egfr.dataset %>%
    mutate(postegfr = ifelse(is.na(postegfr12m), postegfr6m, postegfr12m))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if post eGRF is missing")
    print("################################################")
    print(table(is.na(egfr.dataset$postegfr)))
    print(table(is.na(egfr.dataset$postegfr), egfr.dataset$drugclass))
    
  }
  
  egfr.dataset <- egfr.dataset %>%
    filter(!is.na(postegfr))
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.egfr.dataset <- egfr.dataset %>%
    select(
      # information regarding patient
      patid, pated,
      # response hba1c
      posthba1cfinal,
      # therapies of interest
      drugclass,
      # Sociodemographic features
      agetx, sex, t2dmduration, 
      # Diabetes treatment 
      drugline, ncurrtx, hba1cmonth,
      # Biomarkers
      prehba1c, prebmi, preegfr, preacr, prealbuminblood, prealt, preast, prebilirubin,
      prehaematocrit, prehaemoglobin, prehdl, premap, pretotalcholesterol, pretriglyceride,
      # Comorbidities
      preangina, precld, prediabeticnephropathy, preheartfailure, prehypertension, preihd, premyocardialinfarction, 
      preneuropathy, prepad, preretinopathy, prerevasc, prestroke, pretia, preaf,
      # eGFR analysis
      postegfr
    ) %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### eGFR model - final")
    print("################################################")
    print(nrow(final.egfr.dataset))
    print(table(final.egfr.dataset$drugclass))
    
  }
  
  if (dataset.type == "egfr.dataset") {
    return(final.egfr.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################### CKD outcome population ##############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CKD outcome model")
    print("################################################")
    
  }
  
  ckd.dataset <- final.dataset
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(ckd.dataset$multi_drug_start))
    print(table(ckd.dataset$multi_drug_start, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(ckd.dataset$timeprevcombo_less61))
    print(table(ckd.dataset$timeprevcombo_less61, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(ckd.dataset$hb_extreme_53))
    print(table(ckd.dataset$hb_extreme_53, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(ckd.dataset$prehba1c)))
    print(table(is.na(ckd.dataset$prehba1c), ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(!is.na(prehba1c))
  
  
  ################################################
  ##### Drop if CKD
  ################################################
  
  ckd.dataset <- ckd.dataset %>%
    mutate(no_ckd = ifelse(!is.na(preckd) & (preckd=="stage_3a" | preckd=="stage_3b" | preckd=="stage_4"), 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if CKD")
    print("################################################")
    print(table(ckd.dataset$no_ckd))
    print(table(ckd.dataset$no_ckd, ckd.dataset$drugclass))
    
  }
  
  ckd.dataset <- ckd.dataset %>%
    filter(is.na(no_ckd))
  
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(ckd.dataset$TZD, ckd.dataset$drugclass))
    print("GLP1 treated")
    print(table(ckd.dataset$GLP1, ckd.dataset$drugclass))
    print("SGLT2 treated")
    print(table(ckd.dataset$SGLT2, ckd.dataset$drugclass))
    
  }
  
  
  ckd.dataset <- ckd.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  ckd.dataset <- ckd.dataset %>% 
    mutate(postdrug_stage345 = pmin(postckdstage345date, na.rm = TRUE),
           postdrug_egfr40_or_ckd5 = pmin(egfr40_or_ckd5, na.rm = TRUE)) %>%
    
    
    mutate(postdrug_stage345_censdate = if_else(drugclass=="GLP1",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_sglt2_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_stage345, na.rm=TRUE),
                                                if_else(drugclass=="SGLT2",
                                                        pmin(five_years_post_dstart,
                                                             death_date,
                                                             next_glp1_start,
                                                             next_tzd_start,
                                                             gp_record_end,
                                                             postdrug_stage345, na.rm=TRUE),
                                                        as.Date(NA))),
           
           postdrug_egfr40_or_ckd5_censdate = if_else(drugclass=="GLP1",
                                                      pmin(five_years_post_dstart,
                                                           death_date,
                                                           next_sglt2_start,
                                                           next_tzd_start,
                                                           gp_record_end,
                                                           postdrug_egfr40_or_ckd5, na.rm=TRUE),
                                                      if_else(drugclass=="SGLT2",
                                                              pmin(five_years_post_dstart,
                                                                   death_date,
                                                                   next_glp1_start,
                                                                   next_tzd_start,
                                                                   gp_record_end,
                                                                   postdrug_egfr40_or_ckd5, na.rm=TRUE),
                                                              as.Date(NA))),
           
           
           postdrug_stage345_censvar=ifelse(!is.na(postdrug_stage345) & postdrug_stage345_censdate==postdrug_stage345, 1, 0),
           
           postdrug_egfr40_or_ckd5_censvar=ifelse(!is.na(postdrug_egfr40_or_ckd5) & postdrug_egfr40_or_ckd5_censdate==postdrug_egfr40_or_ckd5, 1, 0),
           
           postdrug_stage345_censtime_yrs=as.numeric(difftime(postdrug_stage345_censdate, dstartdate, unit="days"))/365.25,
           
           postdrug_egfr40_or_ckd5_censtime_yrs=as.numeric(difftime(postdrug_egfr40_or_ckd5_censdate, dstartdate, unit="days"))/365.25
           
    )
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.ckd.dataset <- ckd.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CKD model - final")
    print("################################################")
    print(nrow(final.ckd.dataset))
    print(table(final.ckd.dataset$drugclass))
    
  }
  
  if (dataset.type == "ckd.dataset") {
    return(final.ckd.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ######################### CVD outcome population ##############################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD outcome model")
    print("################################################")
    
  }
  
  cvd.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(cvd.dataset$multi_drug_start))
    print(table(cvd.dataset$multi_drug_start, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  cvd.dataset <- cvd.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(cvd.dataset$timeprevcombo_less61))
    print(table(cvd.dataset$timeprevcombo_less61, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  cvd.dataset <- cvd.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(cvd.dataset$hb_extreme_53))
    print(table(cvd.dataset$hb_extreme_53, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(cvd.dataset$prehba1c)))
    print(table(is.na(cvd.dataset$prehba1c), cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if CVD
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if CVD")
    print("################################################")
    print(table(cvd.dataset$predrug_cvd))
    print(table(cvd.dataset$predrug_cvd, cvd.dataset$drugclass))
    
  }
  
  cvd.dataset <- cvd.dataset %>%
    filter(predrug_cvd == "No")
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(cvd.dataset$TZD, cvd.dataset$drugclass))
    print("GLP1 treated")
    print(table(cvd.dataset$GLP1, cvd.dataset$drugclass))
    print("SGLT2 treated")
    print(table(cvd.dataset$SGLT2, cvd.dataset$drugclass))
    
  }
  
  
  cvd.dataset <- cvd.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  
  ################################################
  ##### Outcome variables
  ################################################
  
  ## MACE: narrow definition (hospitalisation/death): narrow ('incident') MI/stroke HES codes (primary cause only) + CV death in ONS death (primary cause only)
  
  cvd.dataset <- cvd.dataset %>%
    
    mutate(postdrug_mace=pmin(postdrug_first_primary_incident_mi, postdrug_first_primary_incident_stroke, cv_death_date_primary_cause, na.rm=TRUE)) %>%
    
    mutate(postdrug_mace_censdate=if_else(drugclass=="GLP1",
                                          pmin(five_years_post_dstart,
                                               death_date,
                                               next_sglt2_start,
                                               next_tzd_start,
                                               gp_record_end,
                                               postdrug_mace, na.rm=TRUE),
                                          
                                          if_else(drugclass=="SGLT2",
                                                  pmin(five_years_post_dstart,
                                                       death_date,
                                                       next_glp1_start,
                                                       next_tzd_start,
                                                       gp_record_end,
                                                       postdrug_mace, na.rm=TRUE),
                                                  as.Date(NA))),
           
           postdrug_mace_censvar=ifelse(!is.na(postdrug_mace) & postdrug_mace_censdate==postdrug_mace, 1, 0),
           
           postdrug_mace_censtime_yrs=as.numeric(difftime(postdrug_mace_censdate, dstartdate, unit="days"))/365.25)
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.cvd.dataset <- cvd.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.cvd.dataset))
    print(table(final.cvd.dataset$drugclass))
    
  }
  
  if (dataset.type == "cvd.dataset") {
    return(final.cvd.dataset)
  }
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ##################### Heart Failure outcome population ########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Heart Failure outcome model")
    print("################################################")
    
  }
  
  hf.dataset <- final.dataset
  
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(hf.dataset$multi_drug_start))
    print(table(hf.dataset$multi_drug_start, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  hf.dataset <- hf.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(hf.dataset$timeprevcombo_less61))
    print(table(hf.dataset$timeprevcombo_less61, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  hf.dataset <- hf.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(hf.dataset$hb_extreme_53))
    print(table(hf.dataset$hb_extreme_53, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(hf.dataset$prehba1c)))
    print(table(is.na(hf.dataset$prehba1c), hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if heart failure
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if heart failure")
    print("################################################")
    print(table(hf.dataset$preheartfailure))
    print(table(hf.dataset$preheartfailure, hf.dataset$drugclass))
    
  }
  
  hf.dataset <- hf.dataset %>%
    filter(preheartfailure == "No")
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(hf.dataset$TZD, hf.dataset$drugclass))
    print("GLP1 treated")
    print(table(hf.dataset$GLP1, hf.dataset$drugclass))
    print("SGLT2 treated")
    print(table(hf.dataset$SGLT2, hf.dataset$drugclass))
    
  }
  
  
  hf.dataset <- hf.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  hf.dataset <- hf.dataset %>%
    
    ## HF: narrow definition (hospitalisation/death): HF HES codes (primary cause only) + HF death in ONS death (primary cause only)
    
    
    mutate(postdrug_hf=pmin(postdrug_first_primary_hhf, hf_death_date_primary_cause, na.rm=TRUE)) %>%
    
    mutate(postdrug_hf_censdate=if_else(drugclass=="GLP1",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_tzd_start,
                                             gp_record_end,
                                             postdrug_hf, na.rm=TRUE),
                                        
                                        if_else(drugclass=="SGLT2",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_glp1_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_hf, na.rm=TRUE),
                                                as.Date(NA))),
           
           postdrug_hf_censvar=ifelse(!is.na(postdrug_hf) & postdrug_hf_censdate==postdrug_hf, 1, 0),
           
           postdrug_hf_censtime_yrs=as.numeric(difftime(postdrug_hf_censdate, dstartdate, unit="days"))/365.25)
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.hf.dataset <- hf.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.hf.dataset))
    print(table(final.hf.dataset$drugclass))
    
  }
  
  if (dataset.type == "hf.dataset") {
    return(final.hf.dataset)
  }
  
  
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ################# No comorbidities outcome population ########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### No comorbidities outcome model")
    print("################################################")
    
  }
  
  no_co.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(no_co.dataset$multi_drug_start))
    print(table(no_co.dataset$multi_drug_start, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(no_co.dataset$timeprevcombo_less61))
    print(table(no_co.dataset$timeprevcombo_less61, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(no_co.dataset$hb_extreme_53))
    print(table(no_co.dataset$hb_extreme_53, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(no_co.dataset$prehba1c)))
    print(table(is.na(no_co.dataset$prehba1c), no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if comorbidity
  ################################################
  
  no_co.dataset <- no_co.dataset %>%
    mutate(comorbidities = ifelse(predrug_cvd == "No" & preheartfailure == "No" & (is.na(preckd) | (preckd!="stage_3a" & preckd!="stage_3b" & preckd!="stage_4")), NA_real_, 1))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if comorbidities")
    print("################################################")
    print(table(no_co.dataset$comorbidities, no_co.dataset$drugclass))
    
  }
  
  no_co.dataset <- no_co.dataset %>%
    filter(is.na(comorbidities))
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(no_co.dataset$TZD, no_co.dataset$drugclass))
    print("GLP1 treated")
    print(table(no_co.dataset$GLP1, no_co.dataset$drugclass))
    print("SGLT2 treated")
    print(table(no_co.dataset$SGLT2, no_co.dataset$drugclass))
    
  }
  
  
  no_co.dataset <- no_co.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome variables
  ################################################
  
  
  no_co.dataset <- no_co.dataset %>%
    
    
    mutate(postdrug_mace=pmin(postdrug_first_myocardialinfarction, postdrug_first_stroke, cv_death_date_any_cause, na.rm=TRUE),
           postdrug_hf=pmin(postdrug_first_heartfailure, hf_death_date_any_cause, na.rm=TRUE),
           postdrug_stage345 = pmin(postckdstage345date, na.rm = TRUE),
           postdrug_egfr40_or_ckd5 = pmin(egfr40_or_ckd5, na.rm = TRUE)) %>%
    
    
    # Mace: broad definition: MI/stroke GP codes + broad MI/stroke HES codes (any cause) + CV death in ONS death (any cause)
    
    
    mutate(postdrug_mace_censdate=if_else(drugclass=="GLP1",
                                          pmin(five_years_post_dstart,
                                               death_date,
                                               next_sglt2_start,
                                               next_tzd_start,
                                               gp_record_end,
                                               postdrug_mace, na.rm=TRUE),
                                          
                                          if_else(drugclass=="SGLT2",
                                                  pmin(five_years_post_dstart,
                                                       death_date,
                                                       next_glp1_start,
                                                       next_tzd_start,
                                                       gp_record_end,
                                                       postdrug_mace, na.rm=TRUE),
                                                  as.Date(NA))),
           
           postdrug_mace_censvar=ifelse(!is.na(postdrug_mace) & postdrug_mace_censdate==postdrug_mace, 1, 0),
           
           postdrug_mace_censtime_yrs=as.numeric(difftime(postdrug_mace_censdate, dstartdate, unit="days"))/365.25,
           
           
           ## HF: broad definition: HF GP codes + HF HES codes (any cause) + HF death in ONS death (any cause)
           
           
           postdrug_hf_censdate=if_else(drugclass=="GLP1",
                                        pmin(five_years_post_dstart,
                                             death_date,
                                             next_sglt2_start,
                                             next_tzd_start,
                                             gp_record_end,
                                             postdrug_hf, na.rm=TRUE),
                                        
                                        if_else(drugclass=="SGLT2",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_glp1_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_hf, na.rm=TRUE),
                                                as.Date(NA))),
           
           postdrug_hf_censvar=ifelse(!is.na(postdrug_hf) & postdrug_hf_censdate==postdrug_hf, 1, 0),
           
           postdrug_hf_censtime_yrs=as.numeric(difftime(postdrug_hf_censdate, dstartdate, unit="days"))/365.25,
           
           
           
           ## CKD: reaching stage 3/4/5 OR decrease of 40% in egfr and stage 5
           
           
           
           postdrug_stage345_censdate = if_else(drugclass=="GLP1",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_sglt2_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_stage345, na.rm=TRUE),
                                                if_else(drugclass=="SGLT2",
                                                        pmin(five_years_post_dstart,
                                                             death_date,
                                                             next_glp1_start,
                                                             next_tzd_start,
                                                             gp_record_end,
                                                             postdrug_stage345, na.rm=TRUE),
                                                        as.Date(NA))),
           
           postdrug_egfr40_or_ckd5_censdate = if_else(drugclass=="GLP1",
                                                      pmin(five_years_post_dstart,
                                                           death_date,
                                                           next_sglt2_start,
                                                           next_tzd_start,
                                                           gp_record_end,
                                                           postdrug_egfr40_or_ckd5, na.rm=TRUE),
                                                      if_else(drugclass=="SGLT2",
                                                              pmin(five_years_post_dstart,
                                                                   death_date,
                                                                   next_glp1_start,
                                                                   next_tzd_start,
                                                                   gp_record_end,
                                                                   postdrug_egfr40_or_ckd5, na.rm=TRUE),
                                                              as.Date(NA))),
           
           
           postdrug_stage345_censvar=ifelse(!is.na(postdrug_stage345) & postdrug_stage345_censdate==postdrug_stage345, 1, 0),
           
           postdrug_egfr40_or_ckd5_censvar=ifelse(!is.na(postdrug_egfr40_or_ckd5) & postdrug_egfr40_or_ckd5_censdate==postdrug_egfr40_or_ckd5, 1, 0),
           
           postdrug_stage345_censtime_yrs=as.numeric(difftime(postdrug_stage345_censdate, dstartdate, unit="days"))/365.25,
           
           postdrug_egfr40_or_ckd5_censtime_yrs=as.numeric(difftime(postdrug_egfr40_or_ckd5_censdate, dstartdate, unit="days"))/365.25
           
    )
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.no_co.dataset <- no_co.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### CVD model - final")
    print("################################################")
    print(nrow(final.no_co.dataset))
    print(table(final.no_co.dataset$drugclass))
    
  }
  
  if (dataset.type == "no_co.dataset") {
    return(final.no_co.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  ########## No microvascular complications outcome population ##################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### No microvascular complications outcome model")
    print("################################################")
    
  }
  
  micro_comp.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(micro_comp.dataset$multi_drug_start))
    print(table(micro_comp.dataset$multi_drug_start, micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  micro_comp.dataset <- micro_comp.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(micro_comp.dataset$timeprevcombo_less61))
    print(table(micro_comp.dataset$timeprevcombo_less61, micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  micro_comp.dataset <- micro_comp.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(micro_comp.dataset$hb_extreme_53))
    print(table(micro_comp.dataset$hb_extreme_53, micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(micro_comp.dataset$prehba1c)))
    print(table(is.na(micro_comp.dataset$prehba1c), micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(!is.na(prehba1c))
  
  
  
  ################################################
  ##### Drop if comorbidity
  ################################################
  
  micro_comp.dataset <- micro_comp.dataset %>%
    mutate(comorbidities = ifelse(prediabeticnephropathy == "No" & preneuropathy == "No" & preretinopathy == "No", NA_real_, 1))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if comorbidities")
    print("################################################")
    print(table(micro_comp.dataset$comorbidities, micro_comp.dataset$drugclass))
    
  }
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(is.na(comorbidities))
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(micro_comp.dataset$TZD, micro_comp.dataset$drugclass))
    print("GLP1 treated")
    print(table(micro_comp.dataset$GLP1, micro_comp.dataset$drugclass))
    print("SGLT2 treated")
    print(table(micro_comp.dataset$SGLT2, micro_comp.dataset$drugclass))
    
  }
  
  
  micro_comp.dataset <- micro_comp.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  micro_comp.dataset <- micro_comp.dataset %>%
    
    ## Microvascular complications: diabetic nephropathy, retinopathy, neuropahty
    
    
    mutate(postdrug_micro_comp=pmin(postdrug_first_diabeticnephropathy, postdrug_first_neuropathy, postdrug_first_retinopathy, na.rm=TRUE)) %>%
    
    mutate(postdrug_micro_comp_censdate=if_else(drugclass=="GLP1",
                                                pmin(five_years_post_dstart,
                                                     death_date,
                                                     next_sglt2_start,
                                                     next_tzd_start,
                                                     gp_record_end,
                                                     postdrug_micro_comp, na.rm=TRUE),
                                                
                                                if_else(drugclass=="SGLT2",
                                                        pmin(five_years_post_dstart,
                                                             death_date,
                                                             next_glp1_start,
                                                             next_tzd_start,
                                                             gp_record_end,
                                                             postdrug_micro_comp, na.rm=TRUE),
                                                        as.Date(NA))),
           
           postdrug_micro_comp_censvar=ifelse(!is.na(postdrug_micro_comp) & postdrug_micro_comp_censdate==postdrug_micro_comp, 1, 0),
           
           postdrug_micro_comp_censtime_yrs=as.numeric(difftime(postdrug_micro_comp_censdate, dstartdate, unit="days"))/365.25)
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.micro_comp.dataset <- micro_comp.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Microvacular complications model - final")
    print("################################################")
    print(nrow(final.micro_comp.dataset))
    print(table(final.micro_comp.dataset$drugclass))
    
  }
  
  if (dataset.type == "micro_comp.dataset") {
    return(final.micro_comp.dataset)
  }
  
  
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  #################### Retinopathy outcome population ###########################
  ###############################################################################
  #:-----------------------------------------------------------------------------
  ###############################################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Retinopathy outcome model")
    print("################################################")
    
  }
  
  retinopathy.dataset <- final.dataset
  
  
  ################################################
  ##### Drop duplicates (i.e. started treatment on same day)
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop duplicates (i.e. started treatment on same day)")
    print("################################################")
    print(table(retinopathy.dataset$multi_drug_start))
    print(table(retinopathy.dataset$multi_drug_start, retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(multi_drug_start == 0)
  
  
  ################################################
  ##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)
  ################################################
  
  retinopathy.dataset <- retinopathy.dataset %>%
    mutate(timeprevcombo_less61 = ifelse(timeprevcombo <= 61, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if less than 61 days since started previous line of therapy (see Bev finalmerge)")
    print("################################################")
    print(table(retinopathy.dataset$timeprevcombo_less61))
    print(table(retinopathy.dataset$timeprevcombo_less61, retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(is.na(timeprevcombo_less61))
  
  
  ################################################
  ##### Drop if HbA1c <53
  ################################################
  
  retinopathy.dataset <- retinopathy.dataset %>%
    mutate(hb_extreme_53 = ifelse(prehba1c < 53, 1, NA_real_))
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c <53")
    print("################################################")
    print(table(retinopathy.dataset$hb_extreme_53))
    print(table(retinopathy.dataset$hb_extreme_53, retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(is.na(hb_extreme_53))
  
  
  ################################################
  ##### Drop if HbA1c is missing
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if HbA1c is missing")
    print("################################################")
    print(table(is.na(retinopathy.dataset$prehba1c)))
    print(table(is.na(retinopathy.dataset$prehba1c), retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(!is.na(prehba1c))
  
  ################################################
  ##### Drop if retinopathy
  ################################################
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Drop if retinopathy")
    print("################################################")
    print(table(retinopathy.dataset$preretinopathy))
    print(table(retinopathy.dataset$preretinopathy, retinopathy.dataset$drugclass))
    
  }
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(preretinopathy == "No")
  
  
  
  ################################################
  ##### Drop if co-treated with other treatments
  ################################################
  
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("#####  Drop if co-treated with other treatments")
    print("################################################")
    print("TZD treated")
    print(table(retinopathy.dataset$TZD, retinopathy.dataset$drugclass))
    print("GLP1 treated")
    print(table(retinopathy.dataset$GLP1, retinopathy.dataset$drugclass))
    print("SGLT2 treated")
    print(table(retinopathy.dataset$SGLT2, retinopathy.dataset$drugclass))
    
  }
  
  
  retinopathy.dataset <- retinopathy.dataset %>%
    filter(TZD == 0) %>%
    filter(!(GLP1 == 1 & drugclass == "SGLT2")) %>%
    filter(!(SGLT2 == 1 & drugclass == "GLP1"))
  
  
  ################################################
  ##### Outcome
  ################################################
  
  retinopathy.dataset <- retinopathy.dataset %>% 
    mutate(postdrug_retinopathy = pmin(postdrug_first_retinopathy, na.rm = TRUE)) %>%
    
    
    mutate(postdrug_retinopathy_censdate = if_else(drugclass=="GLP1",
                                                   pmin(five_years_post_dstart,
                                                        death_date,
                                                        next_sglt2_start,
                                                        next_tzd_start,
                                                        gp_record_end,
                                                        postdrug_retinopathy, na.rm=TRUE),
                                                   if_else(drugclass=="SGLT2",
                                                           pmin(five_years_post_dstart,
                                                                death_date,
                                                                next_glp1_start,
                                                                next_tzd_start,
                                                                gp_record_end,
                                                                postdrug_retinopathy, na.rm=TRUE),
                                                           as.Date(NA))),
           
           postdrug_retinopathy_censvar=ifelse(!is.na(postdrug_retinopathy) & postdrug_retinopathy_censdate==postdrug_retinopathy, 1, 0),
           
           postdrug_retinopathy_censtime_yrs=as.numeric(difftime(postdrug_retinopathy_censdate, dstartdate, unit="days"))/365.25
           
    )
  
  
  
  #:----------------------------------------------------
  # Select variables needed
  
  final.retinopathy.dataset <- retinopathy.dataset %>%
    as.data.frame()
  
  # printing inclusion patients
  if (dataset.type == "diagnostics") {
    
    print("################################################")
    print("##### Retinopathy outcome model - final")
    print("################################################")
    print(nrow(final.retinopathy.dataset))
    print(table(final.retinopathy.dataset$drugclass))
    
  }
  
  if (dataset.type == "retinopathy.dataset") {
    return(final.retinopathy.dataset)
  }
  
  
  
  
  
}



# This function compares both functions for setting up data and comfirms they are identical
test_functions <- function() {
  
  # Test the synthetic datasets
  synthetic.original <- set_up_data_sglt2_glp1(dataset.type="synthetic")
  synthetic.newer <- set_up_data(dataset.type="synthetic", drugs = c("GLP1", "SGLT2"))
  if(identical(synthetic.original, synthetic.newer) == TRUE) {
    print("The dataset type 'synthetic' is coded correctly")
  } else {
    print("The dataset type 'synthetic' is NOT CODED correctly")
  }
  
  full.cohort.original <- set_up_data_sglt2_glp1(dataset.type="full.cohort")
  full.cohort.newer <- set_up_data(dataset.type="full.cohort", drugs = c("GLP1", "SGLT2"))
  if(identical(full.cohort.original, full.cohort.newer) == TRUE) {
    print("The dataset type 'full.cohort' is coded correctly")
  } else {
    print("The dataset type 'full.cohort' is NOT CODED correctly")
  }
  
  ps.model.train.original <- set_up_data_sglt2_glp1(dataset.type="ps.model.train")
  ps.model.train.newer <- set_up_data(dataset.type="ps.model.train", drugs = c("GLP1", "SGLT2"))
  if(identical(ps.model.train.original, ps.model.train.newer) == TRUE) {
    print("The dataset type 'ps.model.train' is coded correctly")
  } else {
    print("The dataset type 'ps.model.train' is NOT CODED correctly")
  }
  
  ps.model.test.original <- set_up_data_sglt2_glp1(dataset.type="ps.model.test")
  ps.model.test.newer <- set_up_data(dataset.type="ps.model.test", drugs = c("GLP1", "SGLT2"))
  if(identical(ps.model.test.original, ps.model.test.newer) == TRUE) {
    print("The dataset type 'ps.model.test' is coded correctly")
  } else {
    print("The dataset type 'ps.model.test' is NOT CODED correctly")
  }
  
  hba1c.train.original <- set_up_data_sglt2_glp1(dataset.type="hba1c.train")
  hba1c.train.newer <- set_up_data(dataset.type="hba1c.train", drugs = c("GLP1", "SGLT2"))
  if(identical(hba1c.train.original, hba1c.train.newer) == TRUE) {
    print("The dataset type 'hba1c.train' is coded correctly")
  } else {
    print("The dataset type 'hba1c.train' is NOT CODED correctly")
  }
  
  hba1c.test.original <- set_up_data_sglt2_glp1(dataset.type="hba1c.test")
  hba1c.test.newer <- set_up_data(dataset.type="hba1c.test", drugs = c("GLP1", "SGLT2"))
  if(identical(hba1c.test.original, hba1c.test.newer) == TRUE) {
    print("The dataset type 'hba1c.test' is coded correctly")
  } else {
    print("The dataset type 'hba1c.test' is NOT CODED correctly")
  }
  
  weight.dataset.original <- set_up_data_sglt2_glp1(dataset.type="weight.dataset")
  weight.dataset.newer <- set_up_data(dataset.type="weight.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(weight.dataset.original, weight.dataset.newer) == TRUE) {
    print("The dataset type 'weight.dataset' is coded correctly")
  } else {
    print("The dataset type 'weight.dataset' is NOT CODED correctly")
  }
  
  discontinuation.dataset.original <- set_up_data_sglt2_glp1(dataset.type="discontinuation.dataset")
  discontinuation.dataset.newer <- set_up_data(dataset.type="discontinuation.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(discontinuation.dataset.original, discontinuation.dataset.newer) == TRUE) {
    print("The dataset type 'discontinuation.dataset' is coded correctly")
  } else {
    print("The dataset type 'discontinuation.dataset' is NOT CODED correctly")
  }
  
  egfr.dataset.original <- set_up_data_sglt2_glp1(dataset.type="egfr.dataset")
  egfr.dataset.newer <- set_up_data(dataset.type="egfr.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(egfr.dataset.original, egfr.dataset.newer) == TRUE) {
    print("The dataset type 'egfr.dataset' is coded correctly")
  } else {
    print("The dataset type 'egfr.dataset' is NOT CODED correctly")
  }
  
  ckd.dataset.original <- set_up_data_sglt2_glp1(dataset.type="ckd.dataset")
  ckd.dataset.newer <- set_up_data(dataset.type="ckd.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(ckd.dataset.original, ckd.dataset.newer) == TRUE) {
    print("The dataset type 'ckd.dataset' is coded correctly")
  } else {
    print("The dataset type 'ckd.dataset' is NOT CODED correctly")
  }
  
  cvd.dataset.original <- set_up_data_sglt2_glp1(dataset.type="cvd.dataset")
  cvd.dataset.newer <- set_up_data(dataset.type="cvd.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(cvd.dataset.original, cvd.dataset.newer) == TRUE) {
    print("The dataset type 'cvd.dataset' is coded correctly")
  } else {
    print("The dataset type 'cvd.dataset' is NOT CODED correctly")
  }
  
  hf.dataset.original <- set_up_data_sglt2_glp1(dataset.type="hf.dataset")
  hf.dataset.newer <- set_up_data(dataset.type="hf.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(hf.dataset.original, hf.dataset.newer) == TRUE) {
    print("The dataset type 'hf.dataset' is coded correctly")
  } else {
    print("The dataset type 'hf.dataset' is NOT CODED correctly")
  }
  
  no_co.dataset.original <- set_up_data_sglt2_glp1(dataset.type="no_co.dataset")
  no_co.dataset.newer <- set_up_data(dataset.type="no_co.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(no_co.dataset.original, no_co.dataset.newer) == TRUE) {
    print("The dataset type 'no_co.dataset' is coded correctly")
  } else {
    print("The dataset type 'no_co.dataset' is NOT CODED correctly")
  }
  
  micro_comp.dataset.original <- set_up_data_sglt2_glp1(dataset.type="micro_comp.dataset")
  micro_comp.dataset.newer <- set_up_data(dataset.type="micro_comp.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(micro_comp.dataset.original, micro_comp.dataset.newer) == TRUE) {
    print("The dataset type 'micro_comp.dataset' is coded correctly")
  } else {
    print("The dataset type 'micro_comp.dataset' is NOT CODED correctly")
  }
  
  retinopathy.dataset.original <- set_up_data_sglt2_glp1(dataset.type="retinopathy.dataset")
  retinopathy.dataset.newer <- set_up_data(dataset.type="retinopathy.dataset", drugs = c("GLP1", "SGLT2"))
  if(identical(retinopathy.dataset.original, retinopathy.dataset.newer) == TRUE) {
    print("The dataset type 'retinopathy.dataset' is coded correctly")
  } else {
    print("The dataset type 'retinopathy.dataset' is NOT CODED correctly")
  }
  
  
}



