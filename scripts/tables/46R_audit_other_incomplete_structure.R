# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 46R_audit_other_incomplete_structure.R
# Purpose:
#   Identify what remains outside the final frozen cohorts after AKI/covariate
#   exclusions, using final cohort IDs as source of truth.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== AUDIT OTHER INCOMPLETE ANALYTIC STRUCTURE ====")

# -----------------------------------------------------------------------------
# Load upstream and final datasets
# -----------------------------------------------------------------------------

base <- readRDS(file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds")) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id)
  )

aki <- readRDS(file.path(DERIVED_DIR, "11N_aki_outcome_creatinine_7d.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id)) %>%
  dplyr::select(
    op_id,
    baseline_creatinine,
    has_postop_creatinine_7d,
    aki_creat_7d,
    duration_ge_120,
    elective,
    abdominal_department,
    department_clean,
    emop
  )

labs <- readRDS(file.path(DERIVED_DIR, "12R_baseline_haemoglobin_albumin.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id)) %>%
  dplyr::select(op_id, baseline_haemoglobin, baseline_albumin)

flags <- readRDS(file.path(DERIVED_DIR, "24N_operation_level_analysis_flags.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id)) %>%
  dplyr::select(
    op_id,
    any_eligible_interval,
    primary_analysis_op,
    sensitivity_abdominal_op
  )

long_final <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

full_final <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_full_major_noncardiac.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

abd_final <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_elective_abdominal_major.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

# -----------------------------------------------------------------------------
# Final longitudinal availability
# -----------------------------------------------------------------------------

long_ids <- long_final %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    has_longitudinal_rows = dplyr::n() > 0,
    has_final_weight = any(is.finite(final_sw_joint_trunc_1_99)),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# Build audit dataset
# -----------------------------------------------------------------------------

op <- base %>%
  dplyr::left_join(aki, by = "op_id") %>%
  dplyr::left_join(labs, by = "op_id") %>%
  dplyr::left_join(flags, by = "op_id") %>%
  dplyr::left_join(long_ids, by = "op_id")

# -----------------------------------------------------------------------------
# Harmonise duplicated columns after joins
# -----------------------------------------------------------------------------

pick_first_existing <- function(df, candidates) {
  found <- candidates[candidates %in% names(df)]
  if (length(found) == 0) return(NA_character_)
  found[1]
}

duration_col <- pick_first_existing(op, c(
  "duration_ge_120.x",
  "duration_ge_120.y",
  "duration_ge_120"
))

aki_col <- pick_first_existing(op, c(
  "aki_creat_7d.x",
  "aki_creat_7d.y",
  "aki_creat_7d"
))

bmi_col <- pick_first_existing(op, c(
  "bmi.x",
  "bmi.y",
  "bmi",
  "bmi_raw",
  "bmi_derived"
))

asa_col <- pick_first_existing(op, c(
  "asa.x",
  "asa.y",
  "asa"
))

department_col <- pick_first_existing(op, c(
  "department_clean.x",
  "department_clean.y",
  "department_clean",
  "department"
))

# -----------------------------------------------------------------------------
# Create harmonised variables
# -----------------------------------------------------------------------------

if (is.na(duration_col)) {
  op$duration_ge_120_flow <- as.integer(op$intraop_duration_min >= 120)
} else {
  op$duration_ge_120_flow <- as.integer(op[[duration_col]])
}

if (is.na(aki_col)) {
  stop("No AKI outcome column found.")
} else {
  op$aki_creat_7d_flow <- as.integer(op[[aki_col]])
}

if (is.na(bmi_col)) {
  op$bmi_final <- NA_real_
} else {
  op$bmi_final <- op[[bmi_col]]
}

if (is.na(asa_col)) {
  op$asa_final <- NA
} else {
  op$asa_final <- op[[asa_col]]
}

if (is.na(department_col)) {
  op$department_final <- NA_character_
} else {
  op$department_final <- as.character(op[[department_col]])
}

# -----------------------------------------------------------------------------
# Final analytic flags
# -----------------------------------------------------------------------------

op <- op %>%
  dplyr::mutate(
    has_longitudinal_rows = dplyr::coalesce(has_longitudinal_rows, FALSE),
    
    has_final_weight = dplyr::coalesce(has_final_weight, FALSE),
    
    abdominal_candidate = as.integer(
      sensitivity_abdominal_op == 1
    ),
    
    in_full_final = op_id %in% full_final$op_id,
    
    in_abd_final = op_id %in% abd_final$op_id,
    
    aki_evaluable = as.integer(
      !is.na(aki_creat_7d_flow)
    ),
    
    missing_baseline_creatinine = is.na(baseline_creatinine),
    
    missing_baseline_haemoglobin = is.na(baseline_haemoglobin),
    
    missing_baseline_albumin = is.na(baseline_albumin),
    
    missing_bmi = is.na(bmi_final),
    
    missing_asa = is.na(asa_final),
    
    missing_department = is.na(department_final),
    
    incomplete_baseline_covariates =
      missing_baseline_creatinine |
      missing_baseline_haemoglobin |
      missing_baseline_albumin |
      missing_bmi |
      missing_asa |
      missing_department,
    
    unavailable_longitudinal_structure =
      !has_longitudinal_rows |
      !has_final_weight |
      is.na(any_eligible_interval) |
      any_eligible_interval == 0
  )

# -----------------------------------------------------------------------------
# Candidate sets
# -----------------------------------------------------------------------------

full_candidate <- op %>%
  dplyr::filter(duration_ge_120_flow == 1)

abd_candidate <- op %>%
  dplyr::filter(abdominal_candidate == 1)

# -----------------------------------------------------------------------------
# Mutually exclusive reasons
# -----------------------------------------------------------------------------

classify <- function(df, final_col) {
  df %>%
    dplyr::mutate(
      exclusion_reason = dplyr::case_when(
        .data[[final_col]] ~ "Included in final cohort",
        aki_evaluable == 0 ~ "No evaluable postoperative creatinine / AKI outcome",
        missing_baseline_creatinine ~ "Missing baseline creatinine",
        missing_baseline_haemoglobin ~ "Missing baseline haemoglobin",
        missing_baseline_albumin ~ "Missing baseline albumin",
        missing_bmi ~ "Missing BMI",
        missing_asa ~ "Missing ASA physical status",
        missing_department ~ "Missing surgical specialty",
        unavailable_longitudinal_structure ~ "Unavailable eligible longitudinal treatment-weight structure",
        TRUE ~ "Other unresolved discrepancy"
      )
    )
}

breakdown_full <- classify(full_candidate, "in_full_final") %>%
  dplyr::count(exclusion_reason, name = "procedures") %>%
  dplyr::mutate(
    cohort = "Full major non-cardiac candidate cohort",
    denominator = nrow(full_candidate),
    percent = 100 * procedures / denominator
  ) %>%
  dplyr::select(cohort, exclusion_reason, procedures, percent)

breakdown_abd <- classify(abd_candidate, "in_abd_final") %>%
  dplyr::count(exclusion_reason, name = "procedures") %>%
  dplyr::mutate(
    cohort = "Elective major abdominal candidate cohort",
    denominator = nrow(abd_candidate),
    percent = 100 * procedures / denominator
  ) %>%
  dplyr::select(cohort, exclusion_reason, procedures, percent)

breakdown_all <- dplyr::bind_rows(breakdown_full, breakdown_abd) %>%
  dplyr::arrange(cohort, dplyr::desc(procedures))

# -----------------------------------------------------------------------------
# Non-exclusive missingness detail
# -----------------------------------------------------------------------------

detail <- dplyr::bind_rows(
  full_candidate %>%
    dplyr::summarise(
      cohort = "Full major non-cardiac candidate cohort",
      denominator = dplyr::n(),
      no_evaluable_aki = sum(aki_evaluable == 0, na.rm = TRUE),
      missing_baseline_creatinine = sum(missing_baseline_creatinine, na.rm = TRUE),
      missing_baseline_haemoglobin = sum(missing_baseline_haemoglobin, na.rm = TRUE),
      missing_baseline_albumin = sum(missing_baseline_albumin, na.rm = TRUE),
      missing_bmi = sum(missing_bmi, na.rm = TRUE),
      missing_asa = sum(missing_asa, na.rm = TRUE),
      missing_department = sum(missing_department, na.rm = TRUE),
      unavailable_longitudinal_structure = sum(unavailable_longitudinal_structure, na.rm = TRUE)
    ),
  abd_candidate %>%
    dplyr::summarise(
      cohort = "Elective major abdominal candidate cohort",
      denominator = dplyr::n(),
      no_evaluable_aki = sum(aki_evaluable == 0, na.rm = TRUE),
      missing_baseline_creatinine = sum(missing_baseline_creatinine, na.rm = TRUE),
      missing_baseline_haemoglobin = sum(missing_baseline_haemoglobin, na.rm = TRUE),
      missing_baseline_albumin = sum(missing_baseline_albumin, na.rm = TRUE),
      missing_bmi = sum(missing_bmi, na.rm = TRUE),
      missing_asa = sum(missing_asa, na.rm = TRUE),
      missing_department = sum(missing_department, na.rm = TRUE),
      unavailable_longitudinal_structure = sum(unavailable_longitudinal_structure, na.rm = TRUE)
    )
)

# -----------------------------------------------------------------------------
# Key counts
# -----------------------------------------------------------------------------

flow_key <- tibble::tibble(
  metric = c(
    "valid_time_window_procedures",
    "duration_ge_120",
    "full_final",
    "abdominal_candidate_duration_ge_120",
    "abdominal_final"
  ),
  procedures = c(
    dplyr::n_distinct(op$op_id),
    dplyr::n_distinct(op$op_id[op$duration_ge_120_flow == 1]),
    dplyr::n_distinct(full_final$op_id),
    dplyr::n_distinct(op$op_id[op$abdominal_candidate == 1]),
    dplyr::n_distinct(abd_final$op_id)
  )
)

message("==== Flow key ====")
print(as_tibble(flow_key), n = Inf)

message("==== Mutually exclusive exclusion breakdown ====")
print(as_tibble(breakdown_all), n = Inf)

message("==== Non-mutually exclusive missingness detail ====")
print(as_tibble(detail), n = Inf, width = Inf)

readr::write_csv(flow_key, file.path(TABLES_DIR, "46R_flow_key_counts.csv"))
readr::write_csv(breakdown_all, file.path(TABLES_DIR, "46R_flow_exclusion_breakdown_mutually_exclusive.csv"))
readr::write_csv(detail, file.path(TABLES_DIR, "46R_flow_missingness_detail_nonexclusive.csv"))

message("Saved:")
message(file.path(TABLES_DIR, "46R_flow_exclusion_breakdown_mutually_exclusive.csv"))
message("==== AUDIT COMPLETE ====")