# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 45R_flowchart_exclusion_breakdown.R
# Purpose:
#   Decompose flowchart exclusions for:
#     1) full major non-cardiac cohort
#     2) primary elective abdominal cohort
#
#   Produces mutually exclusive exclusion categories:
#     - no evaluable AKI outcome
#     - missing baseline creatinine
#     - missing baseline haemoglobin
#     - missing baseline albumin
#     - missing BMI
#     - missing ASA
#     - missing surgical specialty
#     - missing comorbidity flags
#     - unavailable eligible longitudinal treatment-weight structure
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== FLOWCHART EXCLUSION BREAKDOWN ====")

# -----------------------------------------------------------------------------
# Load upstream and final objects
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
    has_baseline_creatinine,
    has_postop_creatinine_7d,
    aki_creat_7d,
    aki_stage_creat_7d
  )

labs <- readRDS(file.path(DERIVED_DIR, "12R_baseline_haemoglobin_albumin.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id)) %>%
  dplyr::select(
    op_id,
    baseline_haemoglobin,
    baseline_albumin
  )

comorb <- readRDS(file.path(DERIVED_DIR, "14R_baseline_comorbidity_flags.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id)) %>%
  dplyr::select(
    op_id,
    diabetes,
    ckd,
    ischaemic_heart_disease,
    chronic_pulmonary_disease
  )

flags <- readRDS(file.path(DERIVED_DIR, "24N_operation_level_analysis_flags.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

long_final <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

full_final <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_full_major_noncardiac.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

abd_final <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_elective_abdominal_major.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

# -----------------------------------------------------------------------------
# Build upstream operation-level dataset
# -----------------------------------------------------------------------------

op <- base %>%
  dplyr::left_join(aki, by = "op_id") %>%
  dplyr::left_join(labs, by = "op_id") %>%
  dplyr::left_join(comorb, by = "op_id") %>%
  dplyr::left_join(flags, by = "op_id")

# Resolve duplicated columns if present
pick_col <- function(df, candidates) {
  found <- candidates[candidates %in% names(df)]
  if (length(found) == 0) return(NULL)
  found[1]
}

duration_col <- pick_col(op, c("duration_ge_120.x", "duration_ge_120.y", "duration_ge_120"))
aki_col <- pick_col(op, c("aki_creat_7d.x", "aki_creat_7d.y", "aki_creat_7d"))
postop_col <- pick_col(op, c("has_postop_creatinine_7d.x", "has_postop_creatinine_7d.y", "has_postop_creatinine_7d"))
baseline_creat_flag_col <- pick_col(op, c("has_baseline_creatinine.x", "has_baseline_creatinine.y", "has_baseline_creatinine"))
abdominal_col <- pick_col(op, c("abdominal_elective_duration_ge_120", "candidate_abdominal_elective_duration_ge_120", "sensitivity_abdominal_op"))

if (is.null(duration_col)) {
  if ("intraop_duration_min" %in% names(op)) {
    op$duration_ge_120_flow <- as.integer(op$intraop_duration_min >= 120)
  } else if ("operation_duration_min" %in% names(op)) {
    op$duration_ge_120_flow <- as.integer(op$operation_duration_min >= 120)
  } else {
    stop("Cannot define duration_ge_120_flow.")
  }
} else {
  op$duration_ge_120_flow <- as.integer(op[[duration_col]])
}

if (is.null(aki_col)) stop("Cannot identify AKI outcome column.")
op$aki_creat_7d_flow <- as.integer(op[[aki_col]])

if (is.null(postop_col)) {
  op$has_postop_creatinine_7d_flow <- as.integer(!is.na(op$aki_creat_7d_flow))
} else {
  op$has_postop_creatinine_7d_flow <- as.integer(op[[postop_col]])
}

if (is.null(baseline_creat_flag_col)) {
  op$has_baseline_creatinine_flow <- as.integer(!is.na(op$baseline_creatinine))
} else {
  op$has_baseline_creatinine_flow <- as.integer(op[[baseline_creat_flag_col]])
}

if (is.null(abdominal_col)) {
  op$abdominal_elective_duration_ge_120_flow <- as.integer(
    op$duration_ge_120_flow == 1 &
      op$emop == 0 &
      op$department_clean %in% c("GS", "OG", "UR")
  )
} else {
  op$abdominal_elective_duration_ge_120_flow <- as.integer(op[[abdominal_col]])
}

# -----------------------------------------------------------------------------
# Longitudinal eligibility / final weights
# -----------------------------------------------------------------------------

long_op <- long_final %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    has_final_longitudinal_structure = as.integer(dplyr::n() > 0),
    has_final_weight = as.integer(any(is.finite(final_sw_joint_trunc_1_99))),
    .groups = "drop"
  )

op <- op %>%
  dplyr::left_join(long_op, by = "op_id") %>%
  dplyr::mutate(
    has_final_longitudinal_structure = dplyr::coalesce(has_final_longitudinal_structure, 0L),
    has_final_weight = dplyr::coalesce(has_final_weight, 0L),
    
    in_full_final = as.integer(op_id %in% full_final$op_id),
    in_abd_final = as.integer(op_id %in% abd_final$op_id),
    
    aki_evaluable = as.integer(!is.na(aki_creat_7d_flow)),
    
    missing_baseline_creatinine = as.integer(is.na(baseline_creatinine)),
    missing_baseline_haemoglobin = as.integer(is.na(baseline_haemoglobin)),
    missing_baseline_albumin = as.integer(is.na(baseline_albumin)),
    missing_bmi = as.integer(is.na(bmi)),
    missing_asa = as.integer(is.na(asa)),
    missing_surgical_specialty = as.integer(is.na(department_clean)),
    
    missing_diabetes = as.integer(is.na(diabetes)),
    missing_ckd = as.integer(is.na(ckd)),
    missing_ischaemic_heart_disease = as.integer(is.na(ischaemic_heart_disease)),
    missing_chronic_pulmonary_disease = as.integer(is.na(chronic_pulmonary_disease)),
    
    missing_any_comorbidity_flag = as.integer(
      missing_diabetes == 1 |
        missing_ckd == 1 |
        missing_ischaemic_heart_disease == 1 |
        missing_chronic_pulmonary_disease == 1
    ),
    
    missing_any_revised_baseline_covariate = as.integer(
      missing_baseline_creatinine == 1 |
        missing_baseline_haemoglobin == 1 |
        missing_baseline_albumin == 1 |
        missing_bmi == 1 |
        missing_asa == 1 |
        missing_surgical_specialty == 1 |
        missing_any_comorbidity_flag == 1
    ),
    
    unavailable_longitudinal_weight_structure = as.integer(
      has_final_longitudinal_structure == 0 |
        has_final_weight == 0
    )
  )

# -----------------------------------------------------------------------------
# Mutually exclusive exclusion function
# -----------------------------------------------------------------------------

assign_exclusion_reason <- function(df, final_flag_col) {
  df %>%
    dplyr::mutate(
      exclusion_reason = dplyr::case_when(
        .data[[final_flag_col]] == 1 ~ "Included in final cohort",
        
        aki_evaluable == 0 ~ "No evaluable postoperative creatinine / AKI outcome",
        
        missing_baseline_creatinine == 1 ~ "Missing baseline creatinine",
        missing_baseline_haemoglobin == 1 ~ "Missing baseline haemoglobin",
        missing_baseline_albumin == 1 ~ "Missing baseline albumin",
        missing_bmi == 1 ~ "Missing BMI",
        missing_asa == 1 ~ "Missing ASA physical status",
        missing_surgical_specialty == 1 ~ "Missing surgical specialty",
        
        missing_any_comorbidity_flag == 1 ~ "Missing ICD-derived comorbidity flags",
        
        unavailable_longitudinal_weight_structure == 1 ~
          "Unavailable eligible longitudinal treatment-weight structure",
        
        TRUE ~ "Other / unresolved"
      )
    )
}

# -----------------------------------------------------------------------------
# Candidate sets
# -----------------------------------------------------------------------------

full_candidates <- op %>%
  dplyr::filter(duration_ge_120_flow == 1) %>%
  assign_exclusion_reason("in_full_final")

abd_candidates <- op %>%
  dplyr::filter(abdominal_elective_duration_ge_120_flow == 1) %>%
  assign_exclusion_reason("in_abd_final")

# -----------------------------------------------------------------------------
# Breakdown tables
# -----------------------------------------------------------------------------

breakdown_full <- full_candidates %>%
  dplyr::count(exclusion_reason, name = "procedures") %>%
  dplyr::mutate(
    cohort = "Full major non-cardiac candidate cohort",
    denominator = dplyr::n_distinct(full_candidates$op_id),
    percent = 100 * procedures / denominator
  ) %>%
  dplyr::select(cohort, exclusion_reason, procedures, percent)

breakdown_abd <- abd_candidates %>%
  dplyr::count(exclusion_reason, name = "procedures") %>%
  dplyr::mutate(
    cohort = "Elective major abdominal candidate cohort",
    denominator = dplyr::n_distinct(abd_candidates$op_id),
    percent = 100 * procedures / denominator
  ) %>%
  dplyr::select(cohort, exclusion_reason, procedures, percent)

breakdown_all <- dplyr::bind_rows(breakdown_full, breakdown_abd) %>%
  dplyr::arrange(cohort, dplyr::desc(procedures))

# Non-mutually exclusive missingness detail among candidates
missing_detail <- dplyr::bind_rows(
  full_candidates %>%
    dplyr::summarise(
      cohort = "Full major non-cardiac candidate cohort",
      denominator = dplyr::n_distinct(op_id),
      no_evaluable_aki = sum(aki_evaluable == 0, na.rm = TRUE),
      missing_baseline_creatinine = sum(missing_baseline_creatinine == 1, na.rm = TRUE),
      missing_baseline_haemoglobin = sum(missing_baseline_haemoglobin == 1, na.rm = TRUE),
      missing_baseline_albumin = sum(missing_baseline_albumin == 1, na.rm = TRUE),
      missing_bmi = sum(missing_bmi == 1, na.rm = TRUE),
      missing_asa = sum(missing_asa == 1, na.rm = TRUE),
      missing_surgical_specialty = sum(missing_surgical_specialty == 1, na.rm = TRUE),
      missing_any_comorbidity_flag = sum(missing_any_comorbidity_flag == 1, na.rm = TRUE),
      unavailable_longitudinal_weight_structure = sum(unavailable_longitudinal_weight_structure == 1, na.rm = TRUE)
    ),
  abd_candidates %>%
    dplyr::summarise(
      cohort = "Elective major abdominal candidate cohort",
      denominator = dplyr::n_distinct(op_id),
      no_evaluable_aki = sum(aki_evaluable == 0, na.rm = TRUE),
      missing_baseline_creatinine = sum(missing_baseline_creatinine == 1, na.rm = TRUE),
      missing_baseline_haemoglobin = sum(missing_baseline_haemoglobin == 1, na.rm = TRUE),
      missing_baseline_albumin = sum(missing_baseline_albumin == 1, na.rm = TRUE),
      missing_bmi = sum(missing_bmi == 1, na.rm = TRUE),
      missing_asa = sum(missing_asa == 1, na.rm = TRUE),
      missing_surgical_specialty = sum(missing_surgical_specialty == 1, na.rm = TRUE),
      missing_any_comorbidity_flag = sum(missing_any_comorbidity_flag == 1, na.rm = TRUE),
      unavailable_longitudinal_weight_structure = sum(unavailable_longitudinal_weight_structure == 1, na.rm = TRUE)
    )
)

# Flow-compatible counts
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
    dplyr::n_distinct(op$op_id[op$abdominal_elective_duration_ge_120_flow == 1]),
    dplyr::n_distinct(abd_final$op_id)
  )
)

# -----------------------------------------------------------------------------
# Print
# -----------------------------------------------------------------------------

message("==== Flow key ====")
print(flow_key, n = Inf)

message("==== Mutually exclusive exclusion breakdown ====")
print(as_tibble(breakdown_all), n = Inf)

message("==== Non-mutually exclusive missingness detail ====")
print(as_tibble(missing_detail), n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  flow_key,
  file.path(TABLES_DIR, "45R_flow_key_counts.csv")
)

readr::write_csv(
  breakdown_all,
  file.path(TABLES_DIR, "45R_flow_exclusion_breakdown_mutually_exclusive.csv")
)

readr::write_csv(
  missing_detail,
  file.path(TABLES_DIR, "45R_flow_missingness_detail_nonexclusive.csv")
)

saveRDS(
  breakdown_all,
  file.path(DERIVED_DIR, "45R_flow_exclusion_breakdown_mutually_exclusive.rds")
)

saveRDS(
  missing_detail,
  file.path(DERIVED_DIR, "45R_flow_missingness_detail_nonexclusive.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "45R_flow_exclusion_breakdown_mutually_exclusive.csv"))
message(file.path(TABLES_DIR, "45R_flow_missingness_detail_nonexclusive.csv"))
message("==== FLOWCHART EXCLUSION BREAKDOWN COMPLETE ====")