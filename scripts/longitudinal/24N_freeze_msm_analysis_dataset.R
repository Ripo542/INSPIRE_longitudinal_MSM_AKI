# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 24N_freeze_msm_analysis_dataset.R
# Purpose:
#   Merge longitudinal MSM dataset with AKI outcome and freeze primary/sensitivity
#   analysis datasets.
# =============================================================================

SCRIPT <- "24N_freeze_msm_analysis_dataset"

source("R/00_setup/00N_setup_project.R")

message("==== FREEZE MSM ANALYSIS DATASET ====")

# -----------------------------------------------------------------------------
# LOAD
# -----------------------------------------------------------------------------

msm_long <- readRDS(file.path(DERIVED_DIR, "23N_msm_longitudinal_dataset_pre_outcome.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

aki <- readRDS(file.path(DERIVED_DIR, "11N_aki_outcome_creatinine_7d.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

message("MSM long rows: ", nrow(msm_long))
message("MSM long operations: ", dplyr::n_distinct(msm_long$op_id))
message("AKI rows: ", nrow(aki))
message("AKI operations: ", dplyr::n_distinct(aki$op_id))

# -----------------------------------------------------------------------------
# MERGE OUTCOME
# -----------------------------------------------------------------------------

outcome_vars <- aki %>%
  dplyr::select(
    op_id,
    baseline_creatinine,
    baseline_creatinine_source,
    has_baseline_creatinine,
    has_postop_creatinine_7d,
    postop_creatinine_max_7d,
    n_postop_creatinine_7d,
    creatinine_abs_increase_7d,
    creatinine_ratio_7d,
    aki_creat_7d,
    aki_stage_creat_7d
  )

dat <- msm_long %>%
  dplyr::left_join(outcome_vars, by = "op_id")

# -----------------------------------------------------------------------------
# ANALYSIS ELIGIBILITY
# -----------------------------------------------------------------------------

dat <- dat %>%
  dplyr::mutate(
    aki_evaluable = !is.na(aki_creat_7d),
    
    baseline_covariates_complete =
      !is.na(age) &
      !is.na(female) &
      !is.na(asa) &
      !is.na(bmi) &
      !is.na(department_clean) &
      !is.na(baseline_creatinine),
    
    primary_analysis_row =
      primary_duration_ge_120 == TRUE &
      aki_evaluable == TRUE &
      baseline_covariates_complete == TRUE &
      eligible_interval_for_treatment_model == TRUE,
    
    sensitivity_abdominal_row =
      abdominal_elective_duration_ge_120 == TRUE &
      aki_evaluable == TRUE &
      baseline_covariates_complete == TRUE &
      eligible_interval_for_treatment_model == TRUE
  )

primary <- dat %>%
  dplyr::filter(primary_analysis_row)

abdominal <- dat %>%
  dplyr::filter(sensitivity_abdominal_row)

# -----------------------------------------------------------------------------
# OPERATION-LEVEL COHORT QC
# -----------------------------------------------------------------------------

op_flags <- dat %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    duration_ge_120 = dplyr::first(primary_duration_ge_120),
    abdominal_elective_duration_ge_120 = dplyr::first(abdominal_elective_duration_ge_120),
    
    aki_evaluable = dplyr::first(aki_evaluable),
    aki_creat_7d = dplyr::first(aki_creat_7d),
    
    baseline_covariates_complete = dplyr::first(baseline_covariates_complete),
    has_baseline_creatinine = dplyr::first(has_baseline_creatinine),
    has_postop_creatinine_7d = dplyr::first(has_postop_creatinine_7d),
    
    any_eligible_interval = any(eligible_interval_for_treatment_model, na.rm = TRUE),
    
    primary_analysis_op =
      dplyr::first(primary_duration_ge_120) == TRUE &
      dplyr::first(aki_evaluable) == TRUE &
      dplyr::first(baseline_covariates_complete) == TRUE &
      any(eligible_interval_for_treatment_model, na.rm = TRUE),
    
    sensitivity_abdominal_op =
      dplyr::first(abdominal_elective_duration_ge_120) == TRUE &
      dplyr::first(aki_evaluable) == TRUE &
      dplyr::first(baseline_covariates_complete) == TRUE &
      any(eligible_interval_for_treatment_model, na.rm = TRUE),
    
    .groups = "drop"
  )

flow_qc <- tibble::tibble(
  step = c(
    "valid_time_window_non_asa_vi",
    "duration_ge_120",
    "duration_ge_120_aki_evaluable",
    "duration_ge_120_aki_evaluable_baseline_covariates_complete",
    "duration_ge_120_final_primary_with_eligible_intervals",
    "abdominal_elective_duration_ge_120",
    "abdominal_elective_duration_ge_120_aki_evaluable",
    "abdominal_elective_duration_ge_120_final_with_eligible_intervals"
  ),
  n_ops = c(
    dplyr::n_distinct(dat$op_id),
    sum(op_flags$duration_ge_120, na.rm = TRUE),
    sum(op_flags$duration_ge_120 & op_flags$aki_evaluable, na.rm = TRUE),
    sum(op_flags$duration_ge_120 & op_flags$aki_evaluable & op_flags$baseline_covariates_complete, na.rm = TRUE),
    sum(op_flags$primary_analysis_op, na.rm = TRUE),
    sum(op_flags$abdominal_elective_duration_ge_120, na.rm = TRUE),
    sum(op_flags$abdominal_elective_duration_ge_120 & op_flags$aki_evaluable, na.rm = TRUE),
    sum(op_flags$sensitivity_abdominal_op, na.rm = TRUE)
  )
)

message("==== Flow QC ====")
print(flow_qc)

# -----------------------------------------------------------------------------
# FINAL DATASET QC
# -----------------------------------------------------------------------------

primary_qc <- primary %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    aki_events = dplyr::n_distinct(op_id[aki_creat_7d == 1]),
    aki_risk_pct = 100 * aki_events / n_ops,
    
    intervals_median_per_op = stats::median(as.numeric(table(op_id)), na.rm = TRUE),
    intervals_p25_per_op = as.numeric(stats::quantile(as.numeric(table(op_id)), 0.25, na.rm = TRUE)),
    intervals_p75_per_op = as.numeric(stats::quantile(as.numeric(table(op_id)), 0.75, na.rm = TRUE)),
    
    map_missing_pct = 100 * mean(is.na(map)),
    hr_missing_pct = 100 * mean(is.na(hr)),
    
    hypotension65_interval_pct = 100 * mean(hypotension65_current == 1, na.rm = TRUE),
    fluid_interval_pct = 100 * mean(fluid_any == 1, na.rm = TRUE),
    vasopressor_interval_pct = 100 * mean(vasopressor_any == 1, na.rm = TRUE),
    
    ops_any_fluid = dplyr::n_distinct(op_id[fluid_any == 1]),
    ops_any_vasopressor = dplyr::n_distinct(op_id[vasopressor_any == 1]),
    ops_any_eph = dplyr::n_distinct(op_id[eph_any == 1]),
    ops_any_phe = dplyr::n_distinct(op_id[phe_any == 1]),
    ops_any_nepi = dplyr::n_distinct(op_id[nepi_any == 1])
  )

abdominal_qc <- abdominal %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    aki_events = dplyr::n_distinct(op_id[aki_creat_7d == 1]),
    aki_risk_pct = 100 * aki_events / n_ops,
    
    intervals_median_per_op = stats::median(as.numeric(table(op_id)), na.rm = TRUE),
    
    map_missing_pct = 100 * mean(is.na(map)),
    hr_missing_pct = 100 * mean(is.na(hr)),
    
    hypotension65_interval_pct = 100 * mean(hypotension65_current == 1, na.rm = TRUE),
    fluid_interval_pct = 100 * mean(fluid_any == 1, na.rm = TRUE),
    vasopressor_interval_pct = 100 * mean(vasopressor_any == 1, na.rm = TRUE),
    
    ops_any_fluid = dplyr::n_distinct(op_id[fluid_any == 1]),
    ops_any_vasopressor = dplyr::n_distinct(op_id[vasopressor_any == 1])
  )

message("==== Primary MSM dataset QC ====")
print(primary_qc)

message("==== Abdominal sensitivity MSM dataset QC ====")
print(abdominal_qc)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  dat,
  file.path(DERIVED_DIR, "24N_msm_longitudinal_with_outcome_all.rds")
)

saveRDS(
  primary,
  file.path(FROZEN_DIR, "24N_primary_msm_dataset_duration_ge_120.rds")
)

saveRDS(
  abdominal,
  file.path(FROZEN_DIR, "24N_sensitivity_msm_dataset_abdominal_elective_duration_ge_120.rds")
)

saveRDS(
  op_flags,
  file.path(DERIVED_DIR, "24N_operation_level_analysis_flags.rds")
)

readr::write_csv(flow_qc, file.path(QC_DIR, "24N_flow_qc.csv"))
readr::write_csv(primary_qc, file.path(QC_DIR, "24N_primary_msm_dataset_qc.csv"))
readr::write_csv(abdominal_qc, file.path(QC_DIR, "24N_abdominal_msm_dataset_qc.csv"))

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  flow_qc = flow_qc,
  primary_qc = primary_qc,
  abdominal_qc = abdominal_qc
)

write_qc_json(qc, "24N_freeze_msm_analysis_dataset.json")

message("==== MSM ANALYSIS DATASETS FROZEN ====")
message("Primary: ", file.path(FROZEN_DIR, "24N_primary_msm_dataset_duration_ge_120.rds"))
message("Abdominal sensitivity: ", file.path(FROZEN_DIR, "24N_sensitivity_msm_dataset_abdominal_elective_duration_ge_120.rds"))