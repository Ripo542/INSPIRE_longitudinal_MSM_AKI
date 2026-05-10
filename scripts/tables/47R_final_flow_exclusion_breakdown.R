# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 47R_final_flow_exclusion_breakdown.R
# Purpose:
#   Final flowchart counts and exclusion breakdown using the definitive 24R
#   revised longitudinal dataset as source of truth.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== FINAL FLOW EXCLUSION BREAKDOWN ====")

# -----------------------------------------------------------------------------
# Load final revised longitudinal dataset
# -----------------------------------------------------------------------------

long <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id)
  )

base <- readRDS(file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds")) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id)
  )

message("Base valid procedures: ", dplyr::n_distinct(base$op_id))
message("Final revised longitudinal procedures: ", dplyr::n_distinct(long$op_id))

# -----------------------------------------------------------------------------
# Operation-level summary from final 24R dataset
# -----------------------------------------------------------------------------

op <- long %>%
  dplyr::arrange(op_id, interval_id) %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    subject_id = dplyr::first(subject_id),
    
    duration_min = dplyr::first(duration_min),
    full_major_noncardiac = dplyr::first(full_major_noncardiac),
    elective_abdominal_major = dplyr::first(elective_abdominal_major),
    
    revised_primary_full_complete = dplyr::first(revised_primary_full_complete),
    revised_primary_abdominal_complete = dplyr::first(revised_primary_abdominal_complete),
    
    aki_creat_7d = dplyr::first(aki_creat_7d),
    revised_aki_evaluable = dplyr::first(revised_aki_evaluable),
    
    baseline_creatinine = dplyr::first(baseline_creatinine),
    baseline_haemoglobin = dplyr::first(baseline_haemoglobin),
    baseline_albumin = dplyr::first(baseline_albumin),
    bmi = dplyr::first(bmi),
    asa = dplyr::first(asa),
    department_clean = dplyr::first(department_clean),
    
    diabetes = dplyr::first(diabetes),
    ckd = dplyr::first(ckd),
    ischaemic_heart_disease = dplyr::first(ischaemic_heart_disease),
    chronic_pulmonary_disease = dplyr::first(chronic_pulmonary_disease),
    
    n_intervals = dplyr::n(),
    has_final_weight = any(is.finite(final_sw_joint_trunc_1_99)),
    has_lagged_map = any(!is.na(map_lag1)),
    has_lagged_hr = any(!is.na(hr_lag1)),
    has_lagged_fluid = any(!is.na(fluid_any_lag1)),
    has_lagged_vaso = any(!is.na(vasopressor_any_lag1)),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    duration_ge_120 = as.integer(duration_min >= 120),
    
    missing_aki = as.integer(is.na(aki_creat_7d) | revised_aki_evaluable == 0),
    missing_baseline_creatinine = as.integer(is.na(baseline_creatinine)),
    missing_baseline_haemoglobin = as.integer(is.na(baseline_haemoglobin)),
    missing_baseline_albumin = as.integer(is.na(baseline_albumin)),
    missing_bmi = as.integer(is.na(bmi)),
    missing_asa = as.integer(is.na(asa)),
    missing_department = as.integer(is.na(department_clean)),
    
    missing_diabetes = as.integer(is.na(diabetes)),
    missing_ckd = as.integer(is.na(ckd)),
    missing_ihd = as.integer(is.na(ischaemic_heart_disease)),
    missing_pulmonary = as.integer(is.na(chronic_pulmonary_disease)),
    
    unavailable_longitudinal_weight_structure = as.integer(
      n_intervals == 0 |
        !has_final_weight |
        !has_lagged_map |
        !has_lagged_hr |
        !has_lagged_fluid |
        !has_lagged_vaso
    )
  )

# -----------------------------------------------------------------------------
# Helper: mutually exclusive exclusion reason
# -----------------------------------------------------------------------------

classify_reasons <- function(df, final_flag) {
  df %>%
    dplyr::mutate(
      exclusion_reason = dplyr::case_when(
        .data[[final_flag]] == 1 ~ "Included in final cohort",
        
        missing_aki == 1 ~ "No evaluable postoperative creatinine / AKI outcome",
        
        missing_baseline_creatinine == 1 ~ "Missing baseline creatinine",
        missing_baseline_haemoglobin == 1 ~ "Missing baseline haemoglobin",
        missing_baseline_albumin == 1 ~ "Missing baseline albumin",
        missing_bmi == 1 ~ "Missing BMI",
        missing_asa == 1 ~ "Missing ASA physical status",
        missing_department == 1 ~ "Missing surgical specialty",
        
        missing_diabetes == 1 |
          missing_ckd == 1 |
          missing_ihd == 1 |
          missing_pulmonary == 1 ~ "Missing ICD-derived comorbidity flags",
        
        unavailable_longitudinal_weight_structure == 1 ~
          "Unavailable eligible longitudinal treatment-weight structure",
        
        TRUE ~ "Other unresolved discrepancy"
      )
    )
}

# -----------------------------------------------------------------------------
# Candidate cohorts from final 24R flags
# -----------------------------------------------------------------------------

full_candidates <- op %>%
  dplyr::filter(full_major_noncardiac == 1)

abdominal_candidates <- op %>%
  dplyr::filter(elective_abdominal_major == 1)

# -----------------------------------------------------------------------------
# Flow key
# -----------------------------------------------------------------------------

flow_key <- tibble::tibble(
  metric = c(
    "valid_time_window_procedures",
    "duration_ge_120_from_base",
    "final_24r_longitudinal_operations",
    "full_major_noncardiac_candidate",
    "full_major_noncardiac_final",
    "elective_abdominal_candidate",
    "elective_abdominal_final"
  ),
  procedures = c(
    dplyr::n_distinct(base$op_id),
    sum(base$intraop_duration_min >= 120, na.rm = TRUE),
    dplyr::n_distinct(op$op_id),
    dplyr::n_distinct(full_candidates$op_id),
    sum(op$revised_primary_full_complete == 1, na.rm = TRUE),
    dplyr::n_distinct(abdominal_candidates$op_id),
    sum(op$revised_primary_abdominal_complete == 1, na.rm = TRUE)
  ),
  patients = c(
    dplyr::n_distinct(base$subject_id),
    dplyr::n_distinct(base$subject_id[base$intraop_duration_min >= 120]),
    dplyr::n_distinct(op$subject_id),
    dplyr::n_distinct(full_candidates$subject_id),
    dplyr::n_distinct(op$subject_id[op$revised_primary_full_complete == 1]),
    dplyr::n_distinct(abdominal_candidates$subject_id),
    dplyr::n_distinct(op$subject_id[op$revised_primary_abdominal_complete == 1])
  )
)

# -----------------------------------------------------------------------------
# Mutually exclusive breakdown
# -----------------------------------------------------------------------------

breakdown_full <- classify_reasons(
  full_candidates,
  "revised_primary_full_complete"
) %>%
  dplyr::count(exclusion_reason, name = "procedures") %>%
  dplyr::mutate(
    cohort = "Full major non-cardiac candidate cohort",
    denominator = dplyr::n_distinct(full_candidates$op_id),
    percent = 100 * procedures / denominator
  ) %>%
  dplyr::select(cohort, exclusion_reason, procedures, percent)

breakdown_abdominal <- classify_reasons(
  abdominal_candidates,
  "revised_primary_abdominal_complete"
) %>%
  dplyr::count(exclusion_reason, name = "procedures") %>%
  dplyr::mutate(
    cohort = "Elective major abdominal candidate cohort",
    denominator = dplyr::n_distinct(abdominal_candidates$op_id),
    percent = 100 * procedures / denominator
  ) %>%
  dplyr::select(cohort, exclusion_reason, procedures, percent)

breakdown_all <- dplyr::bind_rows(
  breakdown_full,
  breakdown_abdominal
) %>%
  dplyr::arrange(cohort, dplyr::desc(procedures))

# -----------------------------------------------------------------------------
# Non-mutually exclusive detail
# -----------------------------------------------------------------------------

detail <- dplyr::bind_rows(
  full_candidates %>%
    dplyr::summarise(
      cohort = "Full major non-cardiac candidate cohort",
      denominator = dplyr::n_distinct(op_id),
      included_final = sum(revised_primary_full_complete == 1, na.rm = TRUE),
      no_evaluable_aki = sum(missing_aki == 1, na.rm = TRUE),
      missing_baseline_creatinine = sum(missing_baseline_creatinine == 1, na.rm = TRUE),
      missing_baseline_haemoglobin = sum(missing_baseline_haemoglobin == 1, na.rm = TRUE),
      missing_baseline_albumin = sum(missing_baseline_albumin == 1, na.rm = TRUE),
      missing_bmi = sum(missing_bmi == 1, na.rm = TRUE),
      missing_asa = sum(missing_asa == 1, na.rm = TRUE),
      missing_department = sum(missing_department == 1, na.rm = TRUE),
      missing_any_comorbidity_flag = sum(
        missing_diabetes == 1 |
          missing_ckd == 1 |
          missing_ihd == 1 |
          missing_pulmonary == 1,
        na.rm = TRUE
      ),
      unavailable_longitudinal_weight_structure =
        sum(unavailable_longitudinal_weight_structure == 1, na.rm = TRUE)
    ),
  abdominal_candidates %>%
    dplyr::summarise(
      cohort = "Elective major abdominal candidate cohort",
      denominator = dplyr::n_distinct(op_id),
      included_final = sum(revised_primary_abdominal_complete == 1, na.rm = TRUE),
      no_evaluable_aki = sum(missing_aki == 1, na.rm = TRUE),
      missing_baseline_creatinine = sum(missing_baseline_creatinine == 1, na.rm = TRUE),
      missing_baseline_haemoglobin = sum(missing_baseline_haemoglobin == 1, na.rm = TRUE),
      missing_baseline_albumin = sum(missing_baseline_albumin == 1, na.rm = TRUE),
      missing_bmi = sum(missing_bmi == 1, na.rm = TRUE),
      missing_asa = sum(missing_asa == 1, na.rm = TRUE),
      missing_department = sum(missing_department == 1, na.rm = TRUE),
      missing_any_comorbidity_flag = sum(
        missing_diabetes == 1 |
          missing_ckd == 1 |
          missing_ihd == 1 |
          missing_pulmonary == 1,
        na.rm = TRUE
      ),
      unavailable_longitudinal_weight_structure =
        sum(unavailable_longitudinal_weight_structure == 1, na.rm = TRUE)
    )
)

# -----------------------------------------------------------------------------
# Print
# -----------------------------------------------------------------------------

message("==== Flow key ====")
print(tibble::as_tibble(flow_key), n = Inf, width = Inf)

message("==== Mutually exclusive exclusion breakdown ====")
print(tibble::as_tibble(breakdown_all), n = Inf, width = Inf)

message("==== Non-mutually exclusive missingness detail ====")
print(tibble::as_tibble(detail), n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  flow_key,
  file.path(TABLES_DIR, "47R_final_flow_key_counts.csv")
)

readr::write_csv(
  breakdown_all,
  file.path(TABLES_DIR, "47R_final_flow_exclusion_breakdown_mutually_exclusive.csv")
)

readr::write_csv(
  detail,
  file.path(TABLES_DIR, "47R_final_flow_missingness_detail_nonexclusive.csv")
)

saveRDS(
  flow_key,
  file.path(DERIVED_DIR, "47R_final_flow_key_counts.rds")
)

saveRDS(
  breakdown_all,
  file.path(DERIVED_DIR, "47R_final_flow_exclusion_breakdown_mutually_exclusive.rds")
)

saveRDS(
  detail,
  file.path(DERIVED_DIR, "47R_final_flow_missingness_detail_nonexclusive.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "47R_final_flow_key_counts.csv"))
message(file.path(TABLES_DIR, "47R_final_flow_exclusion_breakdown_mutually_exclusive.csv"))
message(file.path(TABLES_DIR, "47R_final_flow_missingness_detail_nonexclusive.csv"))
message("==== FINAL FLOW EXCLUSION BREAKDOWN COMPLETE ====")