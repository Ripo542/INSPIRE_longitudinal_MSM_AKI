# Analytical Framework

## Marginal Structural Models for Time-Varying Confounding

This repository implements marginal structural models with inverse probability weighting to estimate the association between intraoperative fluid and vasopressor administration and postoperative acute kidney injury, adjusting for time-varying confounding affected by prior treatment.

## Study Design

Population: Patients undergoing major abdominal surgery.

Cohort: Defined by inclusion/exclusion criteria documented in scripts/cohort/10N_build_base_cohort.R.

Outcome: Postoperative acute kidney injury defined per KDIGO criteria using serum creatinine. Construction documented in scripts/cohort/11N_build_aki_outcome.R.

Exposure: Fluid administration (crystalloid, colloid) and vasopressor administration (noradrenaline, phenylephrine, others) modeled as time-varying exposures at 5-minute intervals during surgery.

Confounders: Time-fixed baseline characteristics (age, sex, BMI, comorbidities, baseline labs, ASA status) and time-varying intraoperative parameters (vital signs, urine output, cumulative prior treatment).

## Longitudinal Data Structure

The intraoperative period is restructured into 5-minute intervals. For each patient and interval, the following are defined:

- L_t: Lagged confounders (vital signs, fluid balance, prior treatment)
- A_t: Treatment indicators (fluid volume, vasopressor use)
- Y: Binary outcome indicator

Longitudinal restructuring is implemented in scripts/longitudinal/20N_build_longitudinal_intervals.R through scripts/longitudinal/24N_freeze_msm_analysis_dataset.R.

## Inverse Probability Weighting

Treatment allocation models are fit using generalized estimating equations to estimate the probability of observed treatment conditional on lagged confounders (scripts/msm/31N_fit_treatment_models_weights.R). Inverse probability weights are constructed as the ratio of marginal to conditional treatment probabilities.

Weights are accumulated across intervals to the procedure level. Weight truncation is applied and diagnostic assessment performed (scripts/msm/35R_weight_diagnostics.R).

## Outcome Modeling

Outcome Modeling

Weighted logistic regression models were fitted using stabilized inverse probability weights derived from the longitudinal treatment allocation models. Marginal predicted risks and absolute risk differences were estimated from the weighted models across cumulative intraoperative exposure levels.

Model interpretation depends on standard assumptions underlying longitudinal marginal structural models, including exchangeability, positivity, and consistency. As with all observational perioperative analyses, these assumptions cannot be fully verified within the INSPIRE dataset.

## References

Robins JM. Marginal structural models versus structural nested mean models for causal inference. In: Handbook of Lifetable Analysis. Springer, 2001.

Hernán MA, Brumback B, Robins JM. Marginal structural models to estimate the causal effect of zidovudine on the progression of HIV disease. Epidemiology. 2000;11(5):561–570.

Cole SR, Hernán MA. Constructing inverse probability weights for marginal structural models. American Journal of Epidemiology. 2008;168(6):656–664.
