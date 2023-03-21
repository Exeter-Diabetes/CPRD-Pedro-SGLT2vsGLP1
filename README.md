# SGLT2-GLP1
Collection of functions and scripts to investigate how clinical features can be used for prediction.

## Final structure for analysis: (Model 11.5)

1. Fit a BART propensity score model (function _bartMachine::bartMachine_) with all variable available.
2. Perform BART variable selection (function _bartMachine::var_selection_by_permute_) to choose variables for the PS model.
3. Re-fit a BART propensity score model (function _bartMachine::bartMachine_) with the selected variables.
4. Fit a SparseBCF model (function _SparseBCF::SparseBCF_) on complete data and perform variable selection:
    - Only include variable with ~ 20% missingness (discard higher amounts)
    - Select variables with an inclusion proportion above 1/n (n = number of variables used)
5. Fit a BCF model (function _bcf::bcf_) on complete data with the selected variables.
    - Fit two versions of the model:
        - With propensity scores included in the "control" or mu(x) (use include_pi = "control")
        - Without propensity scores included in the model (use include_pi = "none")
    - Compare individual predictions from both models in order to decide on the use of propensity scores.
6. Check model fit for the outcomes:
    - Plot standardised results to check for any structure in the residuals.
7. Check model fit for the treatment effects: (model fitted in observational data)
    - Plot predicted CATE vs ATE for several ntiles of predicted treatment effect:
        - Propensity score matching 1:1 (check whether matched individuals are well balanced)
        - Propensity score matching 1:1 whilst adjusting for all variables used in the BCF model.
        - Adjust for all variables used in the BCF model.



Files:
---
(Developed in CPRD: Aurum download)


- 01: Functions used specifically for this portion.
- 02: Detailed explanation of the selection of cohorts.
- 03: Descriptive analysis of datasets.
- 04: Propensity score model.
- 05: Model heterogeneity.
- 06: Risks/Benefits: hba1c change, eight change, eGFR change, discontinuation, CVD/HF/CKD outcomes, microvascular complications.
- 07: Differential treatment effects.
- 08: Paper plots.
    - .1: Main plots of paper
    - .2: Supplementary plots of paper
    - .3: Plots for DUK.
- 09: Comparison of SGLT2vsGLP1 BCF model to SGLT2vsDPP4 linear model (John Dennis).
- 10: Validation of treatment effects splitting by ethnicity.
- 11: Validation of the excluded individuals that were prescribed semaglutide.
- 12: Validation of treatment effects in those insulin treated.
- 13: Validation of treatment effects in those with/without baseline CVD.


    

