# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 31N_fit_treatment_models_weights.R
# Purpose:
#   Fit treatment models for interval-level fluid and vasopressor exposure,
#   create stabilized inverse probability weights using log-weights, and diagnose
#   operation-level final weights.
# =============================================================================

SCRIPT <- "31N_fit_treatment_models_weights"

source("R/00_setup/00N_setup_project.R")

message("==== FIT TREATMENT MODELS AND CREATE WEIGHTS ====")

# -----------------------------------------------------------------------------
# LOAD
# -----------------------------------------------------------------------------

msm <- readRDS(file.path(DERIVED_DIR, "30N_msm_modeling_dataset.rds"))
formulas <- readRDS(file.path(DERIVED_DIR, "30N_msm_formulas.rds"))

message("Rows loaded: ", nrow(msm))
message("Operations loaded: ", dplyr::n_distinct(msm$op_id))

# -----------------------------------------------------------------------------
# COMPLETE-CASE FOR WEIGHT MODELS
# -----------------------------------------------------------------------------

model_vars <- unique(c(
  all.vars(formulas$fluid),
  all.vars(formulas$vasopressor),
  "op_id", "subject_id", "interval_id", "aki_creat_7d"
))

n_before <- nrow(msm)

msm <- msm %>%
  dplyr::filter(
    dplyr::if_all(
      dplyr::all_of(model_vars),
      ~ !is.na(.x)
    )
  )

message("Rows before complete-case weight filtering: ", n_before)
message("Rows after complete-case weight filtering: ", nrow(msm))
message("Operations after complete-case weight filtering: ", dplyr::n_distinct(msm$op_id))

# -----------------------------------------------------------------------------
# DENOMINATOR MODELS
# -----------------------------------------------------------------------------

message("Fitting denominator model: fluid_any")
fit_fluid_denom <- glm(
  formulas$fluid,
  data = msm,
  family = binomial()
)

message("Fitting denominator model: vasopressor_any")
fit_vaso_denom <- glm(
  formulas$vasopressor,
  data = msm,
  family = binomial()
)

# -----------------------------------------------------------------------------
# NUMERATOR MODELS
# -----------------------------------------------------------------------------
# Stabilized weights: simpler baseline + time models.

formula_fluid_num <- fluid_any ~
  time_elapsed_hr_z +
  age_z +
  female +
  asa +
  bmi_z +
  baseline_creatinine_z +
  department_clean

formula_vaso_num <- vasopressor_any ~
  time_elapsed_hr_z +
  age_z +
  female +
  asa +
  bmi_z +
  baseline_creatinine_z +
  department_clean

message("Fitting numerator model: fluid_any")
fit_fluid_num <- glm(
  formula_fluid_num,
  data = msm,
  family = binomial()
)

message("Fitting numerator model: vasopressor_any")
fit_vaso_num <- glm(
  formula_vaso_num,
  data = msm,
  family = binomial()
)

# -----------------------------------------------------------------------------
# PREDICT PROBABILITIES
# -----------------------------------------------------------------------------

eps <- 1e-6

msm_w <- msm %>%
  dplyr::mutate(
    p_fluid_denom = pmin(pmax(predict(fit_fluid_denom, newdata = msm, type = "response"), eps), 1 - eps),
    p_fluid_num   = pmin(pmax(predict(fit_fluid_num,   newdata = msm, type = "response"), eps), 1 - eps),
    
    p_vaso_denom = pmin(pmax(predict(fit_vaso_denom, newdata = msm, type = "response"), eps), 1 - eps),
    p_vaso_num   = pmin(pmax(predict(fit_vaso_num,   newdata = msm, type = "response"), eps), 1 - eps),
    
    log_sw_fluid_interval = dplyr::if_else(
      fluid_any == 1L,
      log(p_fluid_num) - log(p_fluid_denom),
      log1p(-p_fluid_num) - log1p(-p_fluid_denom)
    ),
    
    log_sw_vaso_interval = dplyr::if_else(
      vasopressor_any == 1L,
      log(p_vaso_num) - log(p_vaso_denom),
      log1p(-p_vaso_num) - log1p(-p_vaso_denom)
    ),
    
    log_sw_joint_interval = log_sw_fluid_interval + log_sw_vaso_interval,
    
    sw_fluid_interval = exp(log_sw_fluid_interval),
    sw_vaso_interval = exp(log_sw_vaso_interval),
    sw_joint_interval = exp(log_sw_joint_interval)
  )

# -----------------------------------------------------------------------------
# CUMULATIVE LOG-WEIGHTS
# -----------------------------------------------------------------------------

msm_w <- msm_w %>%
  dplyr::arrange(op_id, interval_id) %>%
  dplyr::group_by(op_id) %>%
  dplyr::mutate(
    log_sw_fluid_cum = cumsum(log_sw_fluid_interval),
    log_sw_vaso_cum  = cumsum(log_sw_vaso_interval),
    log_sw_joint_cum = cumsum(log_sw_joint_interval)
  ) %>%
  dplyr::ungroup()

# -----------------------------------------------------------------------------
# OPERATION-FINAL WEIGHTS
# -----------------------------------------------------------------------------

op_final_weights <- msm_w %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    final_log_sw_fluid = dplyr::last(log_sw_fluid_cum),
    final_log_sw_vaso = dplyr::last(log_sw_vaso_cum),
    final_log_sw_joint = dplyr::last(log_sw_joint_cum),
    n_intervals = dplyr::n(),
    aki_creat_7d = dplyr::first(aki_creat_7d),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    final_sw_fluid = exp(final_log_sw_fluid),
    final_sw_vaso = exp(final_log_sw_vaso),
    final_sw_joint = exp(final_log_sw_joint)
  )

# -----------------------------------------------------------------------------
# TRUNCATE OPERATION-FINAL JOINT WEIGHTS
# -----------------------------------------------------------------------------

q_1_99 <- stats::quantile(
  op_final_weights$final_sw_joint,
  probs = c(0.01, 0.99),
  na.rm = TRUE
)

q_0_5_99_5 <- stats::quantile(
  op_final_weights$final_sw_joint,
  probs = c(0.005, 0.995),
  na.rm = TRUE
)

op_final_weights <- op_final_weights %>%
  dplyr::mutate(
    final_sw_joint_trunc_1_99 =
      pmin(pmax(final_sw_joint, q_1_99[[1]]), q_1_99[[2]]),
    
    final_sw_joint_trunc_0_5_99_5 =
      pmin(pmax(final_sw_joint, q_0_5_99_5[[1]]), q_0_5_99_5[[2]])
  )

# Attach final operation-level weights back to interval-level dataset
msm_w <- msm_w %>%
  dplyr::left_join(
    op_final_weights %>%
      dplyr::select(
        op_id,
        final_sw_joint,
        final_sw_joint_trunc_1_99,
        final_sw_joint_trunc_0_5_99_5
      ),
    by = "op_id"
  )

# -----------------------------------------------------------------------------
# DIAGNOSTICS
# -----------------------------------------------------------------------------

interval_weight_summary <- msm_w %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    
    sw_fluid_interval_mean = mean(sw_fluid_interval, na.rm = TRUE),
    sw_fluid_interval_p99 = as.numeric(stats::quantile(sw_fluid_interval, 0.99, na.rm = TRUE)),
    sw_fluid_interval_max = max(sw_fluid_interval, na.rm = TRUE),
    
    sw_vaso_interval_mean = mean(sw_vaso_interval, na.rm = TRUE),
    sw_vaso_interval_p99 = as.numeric(stats::quantile(sw_vaso_interval, 0.99, na.rm = TRUE)),
    sw_vaso_interval_max = max(sw_vaso_interval, na.rm = TRUE),
    
    sw_joint_interval_mean = mean(sw_joint_interval, na.rm = TRUE),
    sw_joint_interval_p99 = as.numeric(stats::quantile(sw_joint_interval, 0.99, na.rm = TRUE)),
    sw_joint_interval_max = max(sw_joint_interval, na.rm = TRUE)
  )

op_weight_summary <- op_final_weights %>%
  dplyr::summarise(
    n_ops = dplyr::n(),
    
    final_sw_joint_mean = mean(final_sw_joint, na.rm = TRUE),
    final_sw_joint_p50 = as.numeric(stats::quantile(final_sw_joint, 0.50, na.rm = TRUE)),
    final_sw_joint_p95 = as.numeric(stats::quantile(final_sw_joint, 0.95, na.rm = TRUE)),
    final_sw_joint_p99 = as.numeric(stats::quantile(final_sw_joint, 0.99, na.rm = TRUE)),
    final_sw_joint_p995 = as.numeric(stats::quantile(final_sw_joint, 0.995, na.rm = TRUE)),
    final_sw_joint_max = max(final_sw_joint, na.rm = TRUE),
    
    trunc_1_99_lower = q_1_99[[1]],
    trunc_1_99_upper = q_1_99[[2]],
    final_sw_joint_trunc_1_99_mean = mean(final_sw_joint_trunc_1_99, na.rm = TRUE),
    final_sw_joint_trunc_1_99_p50 = as.numeric(stats::quantile(final_sw_joint_trunc_1_99, 0.50, na.rm = TRUE)),
    final_sw_joint_trunc_1_99_p99 = as.numeric(stats::quantile(final_sw_joint_trunc_1_99, 0.99, na.rm = TRUE)),
    final_sw_joint_trunc_1_99_max = max(final_sw_joint_trunc_1_99, na.rm = TRUE),
    
    trunc_0_5_99_5_lower = q_0_5_99_5[[1]],
    trunc_0_5_99_5_upper = q_0_5_99_5[[2]],
    final_sw_joint_trunc_0_5_99_5_mean = mean(final_sw_joint_trunc_0_5_99_5, na.rm = TRUE),
    final_sw_joint_trunc_0_5_99_5_p995 = as.numeric(stats::quantile(final_sw_joint_trunc_0_5_99_5, 0.995, na.rm = TRUE)),
    final_sw_joint_trunc_0_5_99_5_max = max(final_sw_joint_trunc_0_5_99_5, na.rm = TRUE)
  )

message("==== Interval weight summary ====")
print(interval_weight_summary)

message("==== Operation-final weight summary ====")
print(op_weight_summary)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  msm_w,
  file.path(DERIVED_DIR, "31N_msm_dataset_with_weights.rds")
)

saveRDS(
  op_final_weights,
  file.path(DERIVED_DIR, "31N_operation_final_weights.rds")
)

saveRDS(
  list(
    fluid_denom = fit_fluid_denom,
    fluid_num = fit_fluid_num,
    vaso_denom = fit_vaso_denom,
    vaso_num = fit_vaso_num,
    truncation = list(
      q_1_99 = q_1_99,
      q_0_5_99_5 = q_0_5_99_5
    )
  ),
  file.path(DERIVED_DIR, "31N_treatment_models_and_truncation.rds")
)

readr::write_csv(
  interval_weight_summary,
  file.path(QC_DIR, "31N_interval_weight_summary.csv")
)

readr::write_csv(
  op_weight_summary,
  file.path(QC_DIR, "31N_operation_final_weight_summary.csv")
)

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  interval_weight_summary = interval_weight_summary,
  operation_weight_summary = op_weight_summary
)

write_qc_json(qc, "31N_fit_treatment_models_weights.json")

message("==== TREATMENT MODELS AND WEIGHTS COMPLETE ====")
message("Saved: ", file.path(DERIVED_DIR, "31N_msm_dataset_with_weights.rds"))