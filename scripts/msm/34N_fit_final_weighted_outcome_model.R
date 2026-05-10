# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 34N_fit_final_weighted_outcome_model.R
# Purpose:
#   Final weighted outcome models using survey-weighted quasibinomial inference.
#   Collapses sparse ASA and department levels to avoid leverage instability.
# =============================================================================

SCRIPT <- "34N_fit_final_weighted_outcome_model"

source("R/00_setup/00N_setup_project.R")

message("==== FIT FINAL WEIGHTED OUTCOME MODEL ====")

if (!requireNamespace("survey", quietly = TRUE)) {
  stop("Package 'survey' is required. Install with: install.packages('survey')")
}

suppressPackageStartupMessages({
  library(survey)
})

# -----------------------------------------------------------------------------
# LOAD OPERATION-LEVEL DATASET
# -----------------------------------------------------------------------------

op_dat <- readRDS(file.path(DERIVED_DIR, "32N_operation_level_outcome_dataset.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

message("Operations loaded: ", nrow(op_dat))
message("AKI events: ", sum(op_dat$aki_creat_7d == 1, na.rm = TRUE))
message("AKI risk: ", round(100 * mean(op_dat$aki_creat_7d == 1, na.rm = TRUE), 2), "%")

# -----------------------------------------------------------------------------
# COLLAPSE SPARSE FACTOR LEVELS
# -----------------------------------------------------------------------------

major_departments <- c("GS", "OS", "UR", "CTS", "NS", "OG", "OL", "PS")

op_dat <- op_dat %>%
  dplyr::mutate(
    department_final = dplyr::case_when(
      as.character(department_clean) %in% major_departments ~ as.character(department_clean),
      TRUE ~ "Other"
    ),
    department_final = as.factor(department_final),
    
    asa_final = dplyr::case_when(
      as.character(asa) %in% c("1", "I") ~ "I",
      as.character(asa) %in% c("2", "II") ~ "II",
      as.character(asa) %in% c("3", "III") ~ "III",
      TRUE ~ "IV_V"
    ),
    asa_final = as.factor(asa_final)
  )

# -----------------------------------------------------------------------------
# QC FACTORS
# -----------------------------------------------------------------------------

dept_final_qc <- op_dat %>%
  dplyr::group_by(department_final) %>%
  dplyr::summarise(
    n = dplyr::n(),
    events = sum(aki_creat_7d == 1, na.rm = TRUE),
    risk_pct = 100 * mean(aki_creat_7d == 1, na.rm = TRUE),
    mean_weight_1_99 = mean(sw_1_99, na.rm = TRUE),
    max_weight_1_99 = max(sw_1_99, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(n)

asa_final_qc <- op_dat %>%
  dplyr::group_by(asa_final) %>%
  dplyr::summarise(
    n = dplyr::n(),
    events = sum(aki_creat_7d == 1, na.rm = TRUE),
    risk_pct = 100 * mean(aki_creat_7d == 1, na.rm = TRUE),
    mean_weight_1_99 = mean(sw_1_99, na.rm = TRUE),
    max_weight_1_99 = max(sw_1_99, na.rm = TRUE),
    .groups = "drop"
  )

message("==== Final department QC ====")
print(dept_final_qc)

message("==== Final ASA QC ====")
print(asa_final_qc)

# -----------------------------------------------------------------------------
# FORMULAS
# -----------------------------------------------------------------------------

formula_primary <- aki_creat_7d ~
  hypotension_time_10min +
  map_deficit65_per_10 +
  fluid_intervals_per_5 +
  vasopressor_intervals_per_5 +
  age_z +
  female +
  asa_final +
  bmi_z +
  baseline_creatinine_z +
  department_final +
  n_intervals_z

formula_binary <- aki_creat_7d ~
  hypotension_time_10min +
  map_deficit65_per_10 +
  any_fluid +
  any_vasopressor +
  age_z +
  female +
  asa_final +
  bmi_z +
  baseline_creatinine_z +
  department_final +
  n_intervals_z

# -----------------------------------------------------------------------------
# SURVEY DESIGNS
# -----------------------------------------------------------------------------

des_1_99 <- survey::svydesign(
  ids = ~1,
  weights = ~sw_1_99,
  data = op_dat
)

des_0_5_99_5 <- survey::svydesign(
  ids = ~1,
  weights = ~sw_0_5_99_5,
  data = op_dat
)

des_unweighted <- survey::svydesign(
  ids = ~1,
  weights = ~1,
  data = op_dat
)

# -----------------------------------------------------------------------------
# FIT MODELS
# -----------------------------------------------------------------------------

message("Fitting final primary model: weighted 1–99")
fit_primary_1_99 <- survey::svyglm(
  formula_primary,
  design = des_1_99,
  family = quasibinomial()
)

message("Fitting final sensitivity model: weighted 0.5–99.5")
fit_primary_0_5_99_5 <- survey::svyglm(
  formula_primary,
  design = des_0_5_99_5,
  family = quasibinomial()
)

message("Fitting unweighted comparator")
fit_unweighted <- survey::svyglm(
  formula_primary,
  design = des_unweighted,
  family = quasibinomial()
)

message("Fitting binary-treatment weighted model")
fit_binary_1_99 <- survey::svyglm(
  formula_binary,
  design = des_1_99,
  family = quasibinomial()
)

# -----------------------------------------------------------------------------
# EXTRACT RESULTS
# -----------------------------------------------------------------------------

extract_svy <- function(fit, model_name) {
  sm <- summary(fit)
  ct <- as.data.frame(sm$coefficients)
  ct$term <- rownames(ct)
  
  tibble::tibble(
    model = model_name,
    term = ct$term,
    estimate = ct$Estimate,
    se = ct$`Std. Error`,
    statistic = ct$`t value`,
    p_value = ct$`Pr(>|t|)`,
    or = exp(estimate),
    ci_low = exp(estimate - 1.96 * se),
    ci_high = exp(estimate + 1.96 * se)
  )
}

res_all <- dplyr::bind_rows(
  extract_svy(fit_primary_1_99, "primary_weighted_1_99"),
  extract_svy(fit_primary_0_5_99_5, "sensitivity_weighted_0_5_99_5"),
  extract_svy(fit_unweighted, "unweighted_comparator"),
  extract_svy(fit_binary_1_99, "binary_treatment_weighted_1_99")
)

main_terms <- c(
  "hypotension_time_10min",
  "map_deficit65_per_10",
  "fluid_intervals_per_5",
  "vasopressor_intervals_per_5",
  "any_fluid",
  "any_vasopressor"
)

main_results <- res_all %>%
  dplyr::filter(term %in% main_terms) %>%
  dplyr::mutate(
    estimate_label = sprintf("%.2f (%.2f–%.2f)", or, ci_low, ci_high),
    p_label = dplyr::case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

message("==== Final main results ====")
print(main_results)

# -----------------------------------------------------------------------------
# MODEL STABILITY CHECK
# -----------------------------------------------------------------------------

# Check whether fitted probabilities are reasonable
pred_qc <- tibble::tibble(
  model = c("primary_weighted_1_99", "sensitivity_weighted_0_5_99_5", "unweighted_comparator"),
  pred_min = c(
    min(predict(fit_primary_1_99, type = "response"), na.rm = TRUE),
    min(predict(fit_primary_0_5_99_5, type = "response"), na.rm = TRUE),
    min(predict(fit_unweighted, type = "response"), na.rm = TRUE)
  ),
  pred_p01 = c(
    as.numeric(quantile(predict(fit_primary_1_99, type = "response"), 0.01, na.rm = TRUE)),
    as.numeric(quantile(predict(fit_primary_0_5_99_5, type = "response"), 0.01, na.rm = TRUE)),
    as.numeric(quantile(predict(fit_unweighted, type = "response"), 0.01, na.rm = TRUE))
  ),
  pred_median = c(
    median(predict(fit_primary_1_99, type = "response"), na.rm = TRUE),
    median(predict(fit_primary_0_5_99_5, type = "response"), na.rm = TRUE),
    median(predict(fit_unweighted, type = "response"), na.rm = TRUE)
  ),
  pred_p99 = c(
    as.numeric(quantile(predict(fit_primary_1_99, type = "response"), 0.99, na.rm = TRUE)),
    as.numeric(quantile(predict(fit_primary_0_5_99_5, type = "response"), 0.99, na.rm = TRUE)),
    as.numeric(quantile(predict(fit_unweighted, type = "response"), 0.99, na.rm = TRUE))
  ),
  pred_max = c(
    max(predict(fit_primary_1_99, type = "response"), na.rm = TRUE),
    max(predict(fit_primary_0_5_99_5, type = "response"), na.rm = TRUE),
    max(predict(fit_unweighted, type = "response"), na.rm = TRUE)
  )
)

message("==== Predicted risk QC ====")
print(pred_qc)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  list(
    primary_weighted_1_99 = fit_primary_1_99,
    sensitivity_weighted_0_5_99_5 = fit_primary_0_5_99_5,
    unweighted_comparator = fit_unweighted,
    binary_treatment_weighted_1_99 = fit_binary_1_99
  ),
  file.path(DERIVED_DIR, "34N_final_weighted_outcome_models.rds")
)

saveRDS(
  op_dat,
  file.path(DERIVED_DIR, "34N_final_operation_level_dataset.rds")
)

readr::write_csv(
  res_all,
  file.path(TABLES_DIR, "34N_final_outcome_model_all_terms.csv")
)

readr::write_csv(
  main_results,
  file.path(TABLES_DIR, "34N_final_outcome_model_main_terms.csv")
)

readr::write_csv(
  dept_final_qc,
  file.path(QC_DIR, "34N_department_final_qc.csv")
)

readr::write_csv(
  asa_final_qc,
  file.path(QC_DIR, "34N_asa_final_qc.csv")
)

readr::write_csv(
  pred_qc,
  file.path(QC_DIR, "34N_predicted_risk_qc.csv")
)

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  n_ops = nrow(op_dat),
  aki_events = sum(op_dat$aki_creat_7d == 1, na.rm = TRUE),
  aki_risk_pct = 100 * mean(op_dat$aki_creat_7d == 1, na.rm = TRUE),
  department_final_qc = dept_final_qc,
  asa_final_qc = asa_final_qc,
  predicted_risk_qc = pred_qc,
  main_results = main_results
)

write_qc_json(qc, "34N_fit_final_weighted_outcome_model.json")
# =========================
# SAVE FINAL PRIMARY MODEL
# =========================
saveRDS(
  fit_primary_1_99,
  file = file.path(DERIVED_DIR, "34N_primary_weighted_model.rds")
)
saveRDS(
  list(
    primary_weighted_1_99 = fit_primary_1_99,
    sensitivity_weighted_0_5_99_5 = fit_primary_0_5_99_5,
    unweighted_comparator = fit_unweighted,
    binary_treatment_weighted_1_99 = fit_binary_1_99
  ),
  file = file.path(DERIVED_DIR, "34N_final_weighted_outcome_models.rds")
)
message("==== FINAL WEIGHTED OUTCOME MODEL COMPLETE ====")