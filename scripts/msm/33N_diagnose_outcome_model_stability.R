# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 33N_diagnose_outcome_model_stability.R
# Purpose:
#   Diagnose stability of weighted outcome models:
#   - influential observations
#   - weight distribution among model rows
#   - sparse factor levels
#   - comparison of binomial, quasibinomial, and survey-weighted inference
# =============================================================================

SCRIPT <- "33N_diagnose_outcome_model_stability"

source("R/00_setup/00N_setup_project.R")

message("==== DIAGNOSE OUTCOME MODEL STABILITY ====")

# -----------------------------------------------------------------------------
# REQUIRED PACKAGE
# -----------------------------------------------------------------------------

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
  dplyr::mutate(
    op_id = as.character(op_id),
    asa = as.factor(asa),
    department_clean = as.factor(department_clean)
  )

message("Rows: ", nrow(op_dat))
message("Operations: ", dplyr::n_distinct(op_dat$op_id))
message("AKI events: ", sum(op_dat$aki_creat_7d == 1, na.rm = TRUE))

# -----------------------------------------------------------------------------
# FORMULA
# -----------------------------------------------------------------------------

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

# -----------------------------------------------------------------------------
# SPARSE LEVEL AUDIT
# -----------------------------------------------------------------------------

asa_qc <- op_dat %>%
  dplyr::group_by(asa) %>%
  dplyr::summarise(
    n = dplyr::n(),
    events = sum(aki_creat_7d == 1, na.rm = TRUE),
    risk_pct = 100 * mean(aki_creat_7d == 1, na.rm = TRUE),
    mean_weight = mean(sw_1_99, na.rm = TRUE),
    max_weight = max(sw_1_99, na.rm = TRUE),
    .groups = "drop"
  )

dept_qc <- op_dat %>%
  dplyr::group_by(department_clean) %>%
  dplyr::summarise(
    n = dplyr::n(),
    events = sum(aki_creat_7d == 1, na.rm = TRUE),
    risk_pct = 100 * mean(aki_creat_7d == 1, na.rm = TRUE),
    mean_weight = mean(sw_1_99, na.rm = TRUE),
    max_weight = max(sw_1_99, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(n)

message("==== ASA sparse-level QC ====")
print(asa_qc)

message("==== Department sparse-level QC ====")
print(dept_qc)

# -----------------------------------------------------------------------------
# FIT MODELS FOR COMPARISON
# -----------------------------------------------------------------------------

message("Fitting weighted binomial glm")
fit_binom <- glm(
  formula_primary,
  data = op_dat,
  weights = sw_1_99,
  family = binomial()
)

message("Fitting weighted quasibinomial glm")
fit_quasi <- glm(
  formula_primary,
  data = op_dat,
  weights = sw_1_99,
  family = quasibinomial()
)

message("Fitting survey-weighted quasibinomial model")
des <- survey::svydesign(
  ids = ~1,
  weights = ~sw_1_99,
  data = op_dat
)

fit_svy <- survey::svyglm(
  formula_primary,
  design = des,
  family = quasibinomial()
)

# -----------------------------------------------------------------------------
# INFLUENCE DIAGNOSTICS FROM BINOMIAL MODEL
# -----------------------------------------------------------------------------

hat_vals <- hatvalues(fit_binom)
cook_vals <- cooks.distance(fit_binom)
std_resid <- rstandard(fit_binom, type = "pearson")

diagnostic_dat <- op_dat %>%
  dplyr::mutate(
    row_id = dplyr::row_number(),
    hat = as.numeric(hat_vals),
    cooks_d = as.numeric(cook_vals),
    std_resid = as.numeric(std_resid)
  ) %>%
  dplyr::arrange(dplyr::desc(hat))

top_hat <- diagnostic_dat %>%
  dplyr::select(
    row_id, op_id, aki_creat_7d,
    sw_1_99, sw_0_5_99_5, sw_untruncated,
    hat, cooks_d, std_resid,
    hypotension_time_10min,
    map_deficit65_per_10,
    fluid_intervals_per_5,
    vasopressor_intervals_per_5,
    age, female, asa, bmi,
    baseline_creatinine,
    department_clean,
    n_intervals
  ) %>%
  dplyr::slice_head(n = 20)

top_cook <- diagnostic_dat %>%
  dplyr::arrange(dplyr::desc(cooks_d)) %>%
  dplyr::select(
    row_id, op_id, aki_creat_7d,
    sw_1_99, sw_0_5_99_5, sw_untruncated,
    hat, cooks_d, std_resid,
    hypotension_time_10min,
    map_deficit65_per_10,
    fluid_intervals_per_5,
    vasopressor_intervals_per_5,
    age, female, asa, bmi,
    baseline_creatinine,
    department_clean,
    n_intervals
  ) %>%
  dplyr::slice_head(n = 20)

message("==== Top leverage observations ====")
print(top_hat)

message("==== Top Cook's distance observations ====")
print(top_cook)

# -----------------------------------------------------------------------------
# ROBUST EXTRACTORS
# -----------------------------------------------------------------------------

extract_glm_robust <- function(fit, model_name, vcov_type = "HC3") {
  vc <- sandwich::vcovHC(fit, type = vcov_type)
  ct <- lmtest::coeftest(fit, vcov. = vc)
  
  tibble::tibble(
    model = model_name,
    term = rownames(ct),
    estimate = ct[, "Estimate"],
    se = ct[, "Std. Error"],
    statistic = ct[, "z value"],
    p_value = ct[, "Pr(>|z|)"],
    or = exp(estimate),
    ci_low = exp(estimate - 1.96 * se),
    ci_high = exp(estimate + 1.96 * se)
  )
}

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

res_binom_hc3 <- extract_glm_robust(fit_binom, "weighted_binomial_HC3")
res_quasi_hc3 <- extract_glm_robust(fit_quasi, "weighted_quasibinomial_HC3")
res_svy <- extract_svy(fit_svy, "survey_quasibinomial")

main_terms <- c(
  "hypotension_time_10min",
  "map_deficit65_per_10",
  "fluid_intervals_per_5",
  "vasopressor_intervals_per_5"
)

comparison <- dplyr::bind_rows(
  res_binom_hc3,
  res_quasi_hc3,
  res_svy
) %>%
  dplyr::filter(term %in% main_terms) %>%
  dplyr::mutate(
    estimate_label = sprintf("%.2f (%.2f–%.2f)", or, ci_low, ci_high),
    p_label = dplyr::case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

message("==== Main-term inference comparison ====")
print(comparison)

# -----------------------------------------------------------------------------
# SENSITIVITY: COLLAPSE SPARSE ASA LEVELS IF NEEDED
# -----------------------------------------------------------------------------
# ASA IV/V can be sparse. This version collapses ASA >= IV into a single category.

op_dat_collapsed <- op_dat %>%
  dplyr::mutate(
    asa_collapsed = dplyr::case_when(
      as.character(asa) %in% c("1", "I") ~ "I",
      as.character(asa) %in% c("2", "II") ~ "II",
      as.character(asa) %in% c("3", "III") ~ "III",
      TRUE ~ "IV_or_higher"
    ),
    asa_collapsed = as.factor(asa_collapsed)
  )

formula_collapsed <- aki_creat_7d ~
  hypotension_time_10min +
  map_deficit65_per_10 +
  fluid_intervals_per_5 +
  vasopressor_intervals_per_5 +
  age_z +
  female +
  asa_collapsed +
  bmi_z +
  baseline_creatinine_z +
  department_clean +
  n_intervals_z

fit_quasi_collapsed <- glm(
  formula_collapsed,
  data = op_dat_collapsed,
  weights = sw_1_99,
  family = quasibinomial()
)

res_quasi_collapsed <- extract_glm_robust(
  fit_quasi_collapsed,
  "weighted_quasibinomial_HC3_ASA_collapsed"
) %>%
  dplyr::filter(term %in% main_terms) %>%
  dplyr::mutate(
    estimate_label = sprintf("%.2f (%.2f–%.2f)", or, ci_low, ci_high),
    p_label = dplyr::case_when(
      p_value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p_value)
    )
  )

message("==== ASA-collapsed sensitivity ====")
print(res_quasi_collapsed)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  list(
    fit_binom = fit_binom,
    fit_quasi = fit_quasi,
    fit_svy = fit_svy,
    fit_quasi_collapsed = fit_quasi_collapsed
  ),
  file.path(DERIVED_DIR, "33N_outcome_model_stability_fits.rds")
)

readr::write_csv(asa_qc, file.path(QC_DIR, "33N_asa_sparse_qc.csv"))
readr::write_csv(dept_qc, file.path(QC_DIR, "33N_department_sparse_qc.csv"))
readr::write_csv(top_hat, file.path(QC_DIR, "33N_top_leverage_observations.csv"))
readr::write_csv(top_cook, file.path(QC_DIR, "33N_top_cooks_observations.csv"))
readr::write_csv(comparison, file.path(TABLES_DIR, "33N_main_term_inference_comparison.csv"))
readr::write_csv(res_quasi_collapsed, file.path(TABLES_DIR, "33N_asa_collapsed_sensitivity.csv"))

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  n_ops = nrow(op_dat),
  aki_events = sum(op_dat$aki_creat_7d == 1, na.rm = TRUE),
  asa_qc = asa_qc,
  department_qc = dept_qc,
  comparison = comparison,
  asa_collapsed_sensitivity = res_quasi_collapsed
)

write_qc_json(qc, "33N_diagnose_outcome_model_stability.json")

message("==== OUTCOME MODEL STABILITY DIAGNOSTICS COMPLETE ====")