# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 32N_fit_weighted_outcome_model.R
# Purpose:
#   Fit weighted outcome models for postoperative AKI using truncated stabilized
#   inverse probability weights.
# =============================================================================

SCRIPT <- "32N_fit_weighted_outcome_model"

source("R/00_setup/00N_setup_project.R")

message("==== FIT WEIGHTED OUTCOME MODEL ====")

# -----------------------------------------------------------------------------
# LOAD
# -----------------------------------------------------------------------------

dat <- readRDS(file.path(DERIVED_DIR, "31N_msm_dataset_with_weights.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

message("Rows loaded: ", nrow(dat))
message("Operations loaded: ", dplyr::n_distinct(dat$op_id))

# -----------------------------------------------------------------------------
# COLLAPSE TO OPERATION-LEVEL EXPOSURE SUMMARIES
# -----------------------------------------------------------------------------
# Outcome is operation-level. Treatment/confounder histories were used to build
# weights. The outcome model estimates associations between cumulative treatment/
# hemodynamic exposure and AKI under the weighted pseudo-population.

op_dat <- dat %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    aki_creat_7d = dplyr::first(aki_creat_7d),
    
    # cumulative hemodynamic exposure
    hypotension_intervals = sum(hypotension65_lag1 == 1, na.rm = TRUE),
    map_deficit65_total = max(cumulative_map_deficit65_prev, na.rm = TRUE),
    
    # cumulative treatments
    fluid_intervals = sum(fluid_any == 1, na.rm = TRUE),
    vasopressor_intervals = sum(vasopressor_any == 1, na.rm = TRUE),
    
    any_fluid = as.integer(any(fluid_any == 1, na.rm = TRUE)),
    any_vasopressor = as.integer(any(vasopressor_any == 1, na.rm = TRUE)),
    
    # operation duration / interval count
    n_intervals = dplyr::n(),
    
    # baseline covariates retained for doubly-adjusted outcome model
    age = dplyr::first(age),
    female = dplyr::first(female),
    asa = dplyr::first(asa),
    bmi = dplyr::first(bmi),
    baseline_creatinine = dplyr::first(baseline_creatinine),
    department_clean = dplyr::first(department_clean),
    
    # final weights
    sw_1_99 = dplyr::first(final_sw_joint_trunc_1_99),
    sw_0_5_99_5 = dplyr::first(final_sw_joint_trunc_0_5_99_5),
    sw_untruncated = dplyr::first(final_sw_joint),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    hypotension_time_min = hypotension_intervals * 5,
    fluid_time_min = fluid_intervals * 5,
    vasopressor_time_min = vasopressor_intervals * 5,
    
    # scale interpretable continuous exposures
    hypotension_time_10min = hypotension_time_min / 10,
    map_deficit65_per_10 = map_deficit65_total / 10,
    fluid_intervals_per_5 = fluid_intervals,
    vasopressor_intervals_per_5 = vasopressor_intervals,
    
    age_z = as.numeric(scale(age)),
    bmi_z = as.numeric(scale(bmi)),
    baseline_creatinine_z = as.numeric(scale(baseline_creatinine)),
    n_intervals_z = as.numeric(scale(n_intervals)),
    
    asa = as.factor(asa),
    department_clean = as.factor(department_clean)
  )

message("Operation-level rows: ", nrow(op_dat))
message("AKI events: ", sum(op_dat$aki_creat_7d == 1, na.rm = TRUE))
message("AKI risk: ", round(100 * mean(op_dat$aki_creat_7d == 1, na.rm = TRUE), 2), "%")

# -----------------------------------------------------------------------------
# MODEL FORMULAS
# -----------------------------------------------------------------------------

# Primary weighted MSM-style outcome model.
# We include cumulative hemodynamic burden and cumulative treatment exposures.
# Baseline covariates are included for precision and residual adjustment.

formula_primary <- aki_creat_7d ~
  hypotension_time_10min +
  map_deficit65_per_10 +
  fluid_intervals_per_5 +
  vasopressor_intervals_per_5 +
  age_z +
  female +
  asa +
  bmi_z +
  baseline_creatinine_z +
  department_clean +
  n_intervals_z

formula_binary_treatments <- aki_creat_7d ~
  hypotension_time_10min +
  map_deficit65_per_10 +
  any_fluid +
  any_vasopressor +
  age_z +
  female +
  asa +
  bmi_z +
  baseline_creatinine_z +
  department_clean +
  n_intervals_z

# -----------------------------------------------------------------------------
# FIT MODELS
# -----------------------------------------------------------------------------

message("Fitting primary weighted model: truncation 1–99")
fit_primary_1_99 <- glm(
  formula_primary,
  data = op_dat,
  weights = sw_1_99,
  family = binomial()
)

message("Fitting sensitivity weighted model: truncation 0.5–99.5")
fit_primary_0_5_99_5 <- glm(
  formula_primary,
  data = op_dat,
  weights = sw_0_5_99_5,
  family = binomial()
)

message("Fitting unweighted conventional comparator")
fit_unweighted <- glm(
  formula_primary,
  data = op_dat,
  family = binomial()
)

message("Fitting binary-treatment weighted model")
fit_binary_1_99 <- glm(
  formula_binary_treatments,
  data = op_dat,
  weights = sw_1_99,
  family = binomial()
)

# -----------------------------------------------------------------------------
# ROBUST SE
# -----------------------------------------------------------------------------

robust_or <- function(fit, model_name) {
  vc <- sandwich::vcovHC(fit, type = "HC0")
  ct <- lmtest::coeftest(fit, vcov. = vc)
  
  out <- tibble::tibble(
    model = model_name,
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    robust_se = ct[, "Std. Error"],
    z = ct[, "z value"],
    p_value = ct[, "Pr(>|z|)"],
    or = exp(estimate),
    ci_low = exp(estimate - 1.96 * robust_se),
    ci_high = exp(estimate + 1.96 * robust_se)
  )
  
  out
}

res_primary_1_99 <- robust_or(fit_primary_1_99, "weighted_trunc_1_99")
res_primary_0_5 <- robust_or(fit_primary_0_5_99_5, "weighted_trunc_0_5_99_5")
res_unweighted <- robust_or(fit_unweighted, "unweighted_comparator")
res_binary <- robust_or(fit_binary_1_99, "weighted_binary_treatments_trunc_1_99")

results_all <- dplyr::bind_rows(
  res_primary_1_99,
  res_primary_0_5,
  res_unweighted,
  res_binary
)

# -----------------------------------------------------------------------------
# EXTRACT MAIN TERMS
# -----------------------------------------------------------------------------

main_terms <- c(
  "hypotension_time_10min",
  "map_deficit65_per_10",
  "fluid_intervals_per_5",
  "vasopressor_intervals_per_5",
  "any_fluid",
  "any_vasopressor"
)

main_results <- results_all %>%
  dplyr::filter(term %in% main_terms) %>%
  dplyr::mutate(
    estimate_label = sprintf(
      "%.2f (%.2f–%.2f)",
      or, ci_low, ci_high
    ),
    p_label = dplyr::case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

message("==== Main outcome results ====")
print(main_results)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  list(
    primary_1_99 = fit_primary_1_99,
    primary_0_5_99_5 = fit_primary_0_5_99_5,
    unweighted = fit_unweighted,
    binary_1_99 = fit_binary_1_99
  ),
  file.path(DERIVED_DIR, "32N_weighted_outcome_models.rds")
)

saveRDS(
  op_dat,
  file.path(DERIVED_DIR, "32N_operation_level_outcome_dataset.rds")
)

readr::write_csv(
  results_all,
  file.path(TABLES_DIR, "32N_weighted_outcome_model_all_terms.csv")
)

readr::write_csv(
  main_results,
  file.path(TABLES_DIR, "32N_weighted_outcome_model_main_terms.csv")
)

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  n_ops = nrow(op_dat),
  aki_events = sum(op_dat$aki_creat_7d == 1, na.rm = TRUE),
  aki_risk_pct = 100 * mean(op_dat$aki_creat_7d == 1, na.rm = TRUE),
  main_results = main_results
)

write_qc_json(qc, "32N_fit_weighted_outcome_model.json")

message("==== WEIGHTED OUTCOME MODEL COMPLETE ====")