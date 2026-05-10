# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 44R_flowchart_counts.R
# Purpose:
#   Generate exact counts for revised study flow diagram from upstream cohort.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== FLOWCHART COUNTS FROM UPSTREAM COHORT ====")

# -----------------------------------------------------------------------------
# Load upstream objects
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

message("Base rows: ", nrow(base))
message("Base procedures: ", dplyr::n_distinct(base$op_id))

# -----------------------------------------------------------------------------
# Join
# -----------------------------------------------------------------------------

op <- base %>%
  dplyr::left_join(aki, by = "op_id") %>%
  dplyr::left_join(labs, by = "op_id") %>%
  dplyr::left_join(comorb, by = "op_id") %>%
  dplyr::left_join(flags, by = "op_id")

message("Joined rows: ", nrow(op))
message("Joined procedures: ", dplyr::n_distinct(op$op_id))

# -----------------------------------------------------------------------------
# Helper to choose columns safely
# -----------------------------------------------------------------------------

pick_first_existing <- function(df, candidates) {
  found <- candidates[candidates %in% names(df)]
  if (length(found) == 0) {
    return(NULL)
  }
  found[1]
}

duration_col <- pick_first_existing(
  op,
  c("duration_ge_120.x", "duration_ge_120.y", "duration_ge_120")
)

aki_col <- pick_first_existing(
  op,
  c("aki_creat_7d.x", "aki_creat_7d.y", "aki_creat_7d")
)

postop_creat_col <- pick_first_existing(
  op,
  c("has_postop_creatinine_7d.x", "has_postop_creatinine_7d.y", "has_postop_creatinine_7d")
)

baseline_creat_col <- pick_first_existing(
  op,
  c("has_baseline_creatinine.x", "has_baseline_creatinine.y", "has_baseline_creatinine")
)

abdominal_col <- pick_first_existing(
  op,
  c(
    "abdominal_elective_duration_ge_120",
    "candidate_abdominal_elective_duration_ge_120",
    "sensitivity_abdominal_op"
  )
)

if (is.null(duration_col)) {
  if ("intraop_duration_min" %in% names(op)) {
    op$duration_ge_120_flow <- as.integer(op$intraop_duration_min >= 120)
  } else if ("operation_duration_min" %in% names(op)) {
    op$duration_ge_120_flow <- as.integer(op$operation_duration_min >= 120)
  } else {
    stop("Could not identify duration ≥120 min column.")
  }
} else {
  op$duration_ge_120_flow <- as.integer(op[[duration_col]])
}

if (is.null(aki_col)) {
  stop("Could not identify AKI outcome column.")
} else {
  op$aki_creat_7d_flow <- as.integer(op[[aki_col]])
}

if (is.null(postop_creat_col)) {
  op$has_postop_creatinine_7d_flow <- as.integer(!is.na(op$aki_creat_7d_flow))
} else {
  op$has_postop_creatinine_7d_flow <- as.integer(op[[postop_creat_col]])
}

if (is.null(baseline_creat_col)) {
  op$has_baseline_creatinine_flow <- as.integer(!is.na(op$baseline_creatinine))
} else {
  op$has_baseline_creatinine_flow <- as.integer(op[[baseline_creat_col]])
}

if (is.null(abdominal_col)) {
  if (!all(c("duration_ge_120_flow", "emop", "department_clean") %in% names(op))) {
    stop("Could not define abdominal elective duration ≥120 min flag.")
  }
  
  op$abdominal_elective_duration_ge_120_flow <- as.integer(
    op$duration_ge_120_flow == 1 &
      op$emop == 0 &
      op$department_clean %in% c("GS", "OG", "UR")
  )
} else {
  op$abdominal_elective_duration_ge_120_flow <- as.integer(op[[abdominal_col]])
}

# -----------------------------------------------------------------------------
# Revised complete-case flags
# -----------------------------------------------------------------------------

required_covariates <- c(
  "age",
  "female",
  "asa",
  "bmi",
  "baseline_creatinine",
  "baseline_haemoglobin",
  "baseline_albumin",
  "department_clean",
  "diabetes",
  "ckd",
  "ischaemic_heart_disease",
  "chronic_pulmonary_disease"
)

missing_required <- setdiff(required_covariates, names(op))

if (length(missing_required) > 0) {
  stop(
    "Missing required covariates for flowchart complete-case definition: ",
    paste(missing_required, collapse = ", ")
  )
}

op <- op %>%
  dplyr::mutate(
    aki_evaluable_flow = as.integer(!is.na(aki_creat_7d_flow)),
    
    revised_baseline_covariates_complete = as.integer(
      stats::complete.cases(
        dplyr::across(dplyr::all_of(required_covariates))
      )
    ),
    
    full_revised_complete = as.integer(
      duration_ge_120_flow == 1 &
        aki_evaluable_flow == 1 &
        revised_baseline_covariates_complete == 1
    ),
    
    abdominal_revised_complete = as.integer(
      abdominal_elective_duration_ge_120_flow == 1 &
        aki_evaluable_flow == 1 &
        revised_baseline_covariates_complete == 1
    )
  )

# -----------------------------------------------------------------------------
# Main flow
# -----------------------------------------------------------------------------

flow <- tibble::tibble(
  Step = c(
    "Surgical procedures with valid intraoperative time windows",
    "Procedures lasting ≥120 min",
    "Procedures ≥120 min with evaluable AKI outcome",
    "Full major non-cardiac cohort with complete revised covariates",
    "Elective major abdominal procedures lasting ≥120 min",
    "Primary abdominal cohort with evaluable AKI outcome and complete revised covariates"
  ),
  Procedures = c(
    dplyr::n_distinct(op$op_id),
    dplyr::n_distinct(op$op_id[op$duration_ge_120_flow == 1]),
    dplyr::n_distinct(op$op_id[op$duration_ge_120_flow == 1 & op$aki_evaluable_flow == 1]),
    dplyr::n_distinct(op$op_id[op$full_revised_complete == 1]),
    dplyr::n_distinct(op$op_id[op$abdominal_elective_duration_ge_120_flow == 1]),
    dplyr::n_distinct(op$op_id[op$abdominal_revised_complete == 1])
  ),
  Patients = c(
    dplyr::n_distinct(op$subject_id),
    dplyr::n_distinct(op$subject_id[op$duration_ge_120_flow == 1]),
    dplyr::n_distinct(op$subject_id[op$duration_ge_120_flow == 1 & op$aki_evaluable_flow == 1]),
    dplyr::n_distinct(op$subject_id[op$full_revised_complete == 1]),
    dplyr::n_distinct(op$subject_id[op$abdominal_elective_duration_ge_120_flow == 1]),
    dplyr::n_distinct(op$subject_id[op$abdominal_revised_complete == 1])
  )
)

# -----------------------------------------------------------------------------
# Exclusions
# -----------------------------------------------------------------------------

exclusions <- tibble::tibble(
  Exclusion = c(
    "Duration <120 min",
    "Missing/evaluable AKI unavailable among procedures ≥120 min",
    "Incomplete revised baseline covariates among procedures ≥120 min with evaluable AKI",
    "Not elective major abdominal surgery among procedures ≥120 min",
    "Missing/evaluable AKI or incomplete revised baseline covariates among elective major abdominal procedures ≥120 min",
    "Missing baseline haemoglobin",
    "Missing baseline albumin",
    "Missing baseline creatinine",
    "Missing BMI"
  ),
  Procedures = c(
    sum(op$duration_ge_120_flow == 0, na.rm = TRUE),
    sum(op$duration_ge_120_flow == 1 & op$aki_evaluable_flow == 0, na.rm = TRUE),
    sum(
      op$duration_ge_120_flow == 1 &
        op$aki_evaluable_flow == 1 &
        op$revised_baseline_covariates_complete == 0,
      na.rm = TRUE
    ),
    sum(
      op$duration_ge_120_flow == 1 &
        op$abdominal_elective_duration_ge_120_flow == 0,
      na.rm = TRUE
    ),
    sum(
      op$abdominal_elective_duration_ge_120_flow == 1 &
        op$abdominal_revised_complete == 0,
      na.rm = TRUE
    ),
    sum(is.na(op$baseline_haemoglobin)),
    sum(is.na(op$baseline_albumin)),
    sum(is.na(op$baseline_creatinine)),
    sum(is.na(op$bmi))
  )
)

# -----------------------------------------------------------------------------
# Consistency QC
# -----------------------------------------------------------------------------

expected_full <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_full_major_noncardiac.rds"))
expected_abd <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_elective_abdominal_major.rds"))

qc_consistency <- tibble::tibble(
  Dataset = c(
    "Flow full revised complete",
    "Frozen full operation dataset",
    "Flow abdominal revised complete",
    "Frozen abdominal operation dataset"
  ),
  Procedures = c(
    dplyr::n_distinct(op$op_id[op$full_revised_complete == 1]),
    nrow(expected_full),
    dplyr::n_distinct(op$op_id[op$abdominal_revised_complete == 1]),
    nrow(expected_abd)
  )
)

# -----------------------------------------------------------------------------
# Print and save
# -----------------------------------------------------------------------------

message("==== Flow ====")
print(flow, n = Inf)

message("==== Exclusions ====")
print(exclusions, n = Inf)

message("==== Consistency QC ====")
print(qc_consistency, n = Inf)

readr::write_csv(
  flow,
  file.path(TABLES_DIR, "44R_flowchart_counts.csv")
)

readr::write_csv(
  exclusions,
  file.path(TABLES_DIR, "44R_flowchart_exclusions.csv")
)

readr::write_csv(
  qc_consistency,
  file.path(QC_DIR, "44R_flowchart_consistency_qc.csv")
)

saveRDS(
  flow,
  file.path(DERIVED_DIR, "44R_flowchart_counts.rds")
)

saveRDS(
  exclusions,
  file.path(DERIVED_DIR, "44R_flowchart_exclusions.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "44R_flowchart_counts.csv"))
message(file.path(TABLES_DIR, "44R_flowchart_exclusions.csv"))
message(file.path(QC_DIR, "44R_flowchart_consistency_qc.csv"))
message("==== FLOWCHART COUNTS COMPLETE ====")