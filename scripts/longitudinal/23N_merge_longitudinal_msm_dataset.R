# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 23N_merge_longitudinal_msm_dataset.R
# Purpose:
#   Merge 5-min hemodynamic intervals with 5-min treatment exposures and create
#   time-varying history variables for MSM modelling.
# =============================================================================

SCRIPT <- "23N_merge_longitudinal_msm_dataset"

source("R/00_setup/00N_setup_project.R")

message("==== MERGE LONGITUDINAL MSM DATASET ====")

# -----------------------------------------------------------------------------
# LOAD
# -----------------------------------------------------------------------------

long_hemo <- readRDS(file.path(DERIVED_DIR, "20N_longitudinal_base_intervals.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

tx_long <- readRDS(file.path(DERIVED_DIR, "22N_treatments_5min_long.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

base <- readRDS(file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds")) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id)
  )

message("Hemo rows: ", nrow(long_hemo))
message("Treatment rows: ", nrow(tx_long))
message("Base operations: ", dplyr::n_distinct(base$op_id))

# -----------------------------------------------------------------------------
# MERGE HEMO + TREATMENT
# -----------------------------------------------------------------------------

msm_long <- long_hemo %>%
  dplyr::left_join(
    tx_long %>%
      dplyr::select(
        op_id, t,
        fluid_any, crystalloid_any, colloid_any, blood_any,
        fluid_total_value, crystalloid_total_value, colloid_total_value, blood_total_value,
        vasopressor_any, eph_any, phe_any, nepi_any,
        eph_value, phe_value, nepi_value,
        other_vasoactive_any,
        n_treatment_records
      ),
    by = c("op_id", "t")
  )

# Missing treatments are true non-treatment intervals
msm_long <- msm_long %>%
  dplyr::mutate(
    dplyr::across(
      c(
        fluid_any, crystalloid_any, colloid_any, blood_any,
        vasopressor_any, eph_any, phe_any, nepi_any,
        other_vasoactive_any
      ),
      ~ dplyr::coalesce(.x, 0L)
    ),
    dplyr::across(
      c(
        fluid_total_value, crystalloid_total_value, colloid_total_value,
        colloid_total_value, blood_total_value,
        eph_value, phe_value, nepi_value,
        n_treatment_records
      ),
      ~ dplyr::coalesce(.x, 0)
    )
  )

# -----------------------------------------------------------------------------
# ADD BASELINE / PROCEDURAL VARIABLES
# -----------------------------------------------------------------------------

baseline_vars <- base %>%
  dplyr::select(
    op_id, subject_id,
    age, female, asa, bmi, department_clean, emop,
    t0, t1, intraop_duration_min,
    operation_duration_min, anesthesia_duration_min, or_duration_min,
    abdominal_department, elective, duration_ge_120
  )

msm_long <- msm_long %>%
  dplyr::left_join(baseline_vars, by = c("op_id", "subject_id"))

# -----------------------------------------------------------------------------
# CREATE TIME-VARYING HISTORY VARIABLES
# -----------------------------------------------------------------------------

msm_long <- msm_long %>%
  dplyr::arrange(op_id, t) %>%
  dplyr::group_by(op_id) %>%
  dplyr::mutate(
    interval_id = dplyr::row_number(),
    time_elapsed_min = t,
    time_elapsed_hr = time_elapsed_min / 60,
    time_remaining_min = intraop_duration_min - time_elapsed_min,
    
    # Current hemodynamic state
    hypotension65_current = dplyr::if_else(!is.na(map) & map < 65, 1L, 0L, missing = NA_integer_),
    map_deficit65_current = dplyr::if_else(!is.na(map), pmax(65 - map, 0), NA_real_),
    
    # Lagged physiological state: this is key for treatment models
    map_lag1 = dplyr::lag(map),
    hr_lag1 = dplyr::lag(hr),
    hypotension65_lag1 = dplyr::lag(hypotension65_current),
    map_deficit65_lag1 = dplyr::lag(map_deficit65_current),
    
    # Previous treatment
    fluid_any_lag1 = dplyr::lag(fluid_any, default = 0L),
    vasopressor_any_lag1 = dplyr::lag(vasopressor_any, default = 0L),
    
    # MAP trend based on prior observed state
    map_change_lag1 = map - dplyr::lag(map),
    map_change_prev = dplyr::lag(map_change_lag1),
    
    # Cumulative prior physiologic burden
    cumulative_hypotension_intervals_prev =
      dplyr::lag(cumsum(dplyr::coalesce(hypotension65_current, 0L)), default = 0L),
    
    cumulative_map_deficit65_prev =
      dplyr::lag(cumsum(dplyr::coalesce(map_deficit65_current, 0)), default = 0),
    
    # Cumulative prior treatment exposure
    cumulative_fluid_any_intervals_prev =
      dplyr::lag(cumsum(fluid_any), default = 0L),
    
    cumulative_fluid_value_prev =
      dplyr::lag(cumsum(fluid_total_value), default = 0),
    
    cumulative_vasopressor_intervals_prev =
      dplyr::lag(cumsum(vasopressor_any), default = 0L),
    
    cumulative_eph_intervals_prev =
      dplyr::lag(cumsum(eph_any), default = 0L),
    
    cumulative_phe_intervals_prev =
      dplyr::lag(cumsum(phe_any), default = 0L),
    
    cumulative_nepi_intervals_prev =
      dplyr::lag(cumsum(nepi_any), default = 0L),
    
    # Previous observed coverage
    prior_map_observed_intervals =
      dplyr::lag(cumsum(!is.na(map)), default = 0L),
    
    prior_hr_observed_intervals =
      dplyr::lag(cumsum(!is.na(hr)), default = 0L),
    
    map_observed_fraction_prev =
      prior_map_observed_intervals / pmax(interval_id - 1, 1),
    
    hr_observed_fraction_prev =
      prior_hr_observed_intervals / pmax(interval_id - 1, 1)
  ) %>%
  dplyr::ungroup()

# -----------------------------------------------------------------------------
# DEFINE MODELLING FLAGS
# -----------------------------------------------------------------------------

msm_long <- msm_long %>%
  dplyr::mutate(
    # Analysis starts at interval 2 because lagged physiology is required
    eligible_interval_for_treatment_model =
      interval_id > 1 &
      !is.na(map_lag1) &
      !is.na(hr_lag1) &
      !is.na(hypotension65_lag1),
    
    # Primary broad cohort candidate
    primary_duration_ge_120 =
      duration_ge_120,
    
    # More homogeneous sensitivity cohort
    abdominal_elective_duration_ge_120 =
      abdominal_department == TRUE &
      duration_ge_120 == TRUE &
      (elective == TRUE | is.na(elective))
  )

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

qc_main <- msm_long %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    
    map_missing_pct = 100 * mean(is.na(map)),
    hr_missing_pct = 100 * mean(is.na(hr)),
    
    eligible_treatment_model_rows = sum(eligible_interval_for_treatment_model, na.rm = TRUE),
    eligible_treatment_model_ops = dplyr::n_distinct(op_id[eligible_interval_for_treatment_model]),
    
    primary_duration_ge_120_ops = dplyr::n_distinct(op_id[primary_duration_ge_120]),
    abdominal_elective_duration_ge_120_ops =
      dplyr::n_distinct(op_id[abdominal_elective_duration_ge_120]),
    
    fluid_interval_pct = 100 * mean(fluid_any == 1, na.rm = TRUE),
    vasopressor_interval_pct = 100 * mean(vasopressor_any == 1, na.rm = TRUE),
    
    hypotension65_interval_pct = 100 * mean(hypotension65_current == 1, na.rm = TRUE)
  )

message("==== MSM longitudinal dataset QC ====")
print(qc_main)

op_qc <- msm_long %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    n_intervals = dplyr::n(),
    duration_min = dplyr::first(intraop_duration_min),
    primary_duration_ge_120 = dplyr::first(primary_duration_ge_120),
    abdominal_elective_duration_ge_120 = dplyr::first(abdominal_elective_duration_ge_120),
    
    map_nonmissing_frac = mean(!is.na(map)),
    hr_nonmissing_frac = mean(!is.na(hr)),
    
    any_hypotension65 = any(hypotension65_current == 1, na.rm = TRUE),
    any_fluid = any(fluid_any == 1, na.rm = TRUE),
    any_vasopressor = any(vasopressor_any == 1, na.rm = TRUE),
    any_eph = any(eph_any == 1, na.rm = TRUE),
    any_phe = any(phe_any == 1, na.rm = TRUE),
    any_nepi = any(nepi_any == 1, na.rm = TRUE),
    
    total_fluid_value = sum(fluid_total_value, na.rm = TRUE),
    total_vasopressor_intervals = sum(vasopressor_any, na.rm = TRUE),
    cumulative_map_deficit65_final = max(cumulative_map_deficit65_prev, na.rm = TRUE),
    
    .groups = "drop"
  )

op_summary <- op_qc %>%
  dplyr::summarise(
    n_ops = dplyr::n(),
    n_ops_primary_duration_ge_120 = sum(primary_duration_ge_120, na.rm = TRUE),
    n_ops_abdominal_elective_duration_ge_120 = sum(abdominal_elective_duration_ge_120, na.rm = TRUE),
    
    intervals_median = stats::median(n_intervals, na.rm = TRUE),
    intervals_p25 = as.numeric(stats::quantile(n_intervals, 0.25, na.rm = TRUE)),
    intervals_p75 = as.numeric(stats::quantile(n_intervals, 0.75, na.rm = TRUE)),
    
    map_nonmissing_frac_median = stats::median(map_nonmissing_frac, na.rm = TRUE),
    hr_nonmissing_frac_median = stats::median(hr_nonmissing_frac, na.rm = TRUE),
    
    ops_any_hypotension65 = sum(any_hypotension65, na.rm = TRUE),
    ops_any_fluid = sum(any_fluid, na.rm = TRUE),
    ops_any_vasopressor = sum(any_vasopressor, na.rm = TRUE),
    ops_any_eph = sum(any_eph, na.rm = TRUE),
    ops_any_phe = sum(any_phe, na.rm = TRUE),
    ops_any_nepi = sum(any_nepi, na.rm = TRUE)
  )

message("==== Operation-level MSM QC ====")
print(op_summary)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  msm_long,
  file.path(DERIVED_DIR, "23N_msm_longitudinal_dataset_pre_outcome.rds")
)

saveRDS(
  op_qc,
  file.path(DERIVED_DIR, "23N_msm_operation_qc_pre_outcome.rds")
)

readr::write_csv(
  qc_main,
  file.path(QC_DIR, "23N_msm_longitudinal_qc.csv")
)

readr::write_csv(
  op_summary,
  file.path(QC_DIR, "23N_msm_operation_qc.csv")
)

qc_json <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  longitudinal_qc = qc_main,
  operation_qc = op_summary
)

write_qc_json(qc_json, "23N_merge_longitudinal_msm_dataset.json")

message("==== MSM LONGITUDINAL DATASET COMPLETE ====")
message("Saved: ", file.path(DERIVED_DIR, "23N_msm_longitudinal_dataset_pre_outcome.rds"))