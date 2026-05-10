# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 48R_audit_longitudinal_dataset_losses.R
# Purpose:
#   Decompose exclusions between duration >=120 min cohort and the revised
#   longitudinal analytic dataset (24R).
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== AUDIT LONGITUDINAL DATASET LOSSES ====")

# -----------------------------------------------------------------------------
# Load data
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
    has_baseline_creatinine,
    has_postop_creatinine_7d,
    aki_creat_7d,
    duration_ge_120
  )

flags <- readRDS(file.path(DERIVED_DIR, "24N_operation_level_analysis_flags.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

long24r <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

# -----------------------------------------------------------------------------
# Final longitudinal operation IDs
# -----------------------------------------------------------------------------

long_ids <- long24r %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    n_intervals = dplyr::n(),
    has_weight = any(is.finite(final_sw_joint_trunc_1_99)),
    has_map_lag = any(!is.na(map_lag1)),
    has_hr_lag = any(!is.na(hr_lag1)),
    has_fluid_lag = any(!is.na(fluid_any_lag1)),
    has_vaso_lag = any(!is.na(vasopressor_any_lag1)),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# Build upstream duration >=120 dataset
# -----------------------------------------------------------------------------

op <- base %>%
  dplyr::left_join(aki, by = "op_id") %>%
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

duration_col <- pick_first_existing(
  op,
  c("duration_ge_120.x", "duration_ge_120.y", "duration_ge_120")
)

postop_col <- pick_first_existing(
  op,
  c("has_postop_creatinine_7d.x", "has_postop_creatinine_7d.y", "has_postop_creatinine_7d")
)

baseline_col <- pick_first_existing(
  op,
  c("has_baseline_creatinine.x", "has_baseline_creatinine.y", "has_baseline_creatinine")
)

if (is.na(duration_col)) {
  if (!"intraop_duration_min" %in% names(op)) {
    stop("No duration_ge_120 or intraop_duration_min column found.")
  }
  op$duration_ge_120_flow <- as.integer(op$intraop_duration_min >= 120)
} else {
  op$duration_ge_120_flow <- as.integer(op[[duration_col]])
}

if (is.na(postop_col)) {
  stop("No has_postop_creatinine_7d column found.")
} else {
  op$has_postop_creatinine_7d_flow <- as.integer(op[[postop_col]])
}

if (is.na(baseline_col)) {
  stop("No has_baseline_creatinine column found.")
} else {
  op$has_baseline_creatinine_flow <- as.integer(op[[baseline_col]])
}

# -----------------------------------------------------------------------------
# Restrict to duration >=120
# -----------------------------------------------------------------------------

candidate <- op %>%
  dplyr::filter(duration_ge_120_flow == 1)

message("Duration >=120 candidate operations: ", nrow(candidate))

# -----------------------------------------------------------------------------
# Final longitudinal inclusion and exclusion flags
# -----------------------------------------------------------------------------

candidate <- candidate %>%
  dplyr::mutate(
    included_long24r = op_id %in% long24r$op_id,
    
    missing_postop_creatinine = (
      is.na(has_postop_creatinine_7d_flow) |
        has_postop_creatinine_7d_flow == 0
    ),
    
    missing_baseline_creatinine = (
      is.na(has_baseline_creatinine_flow) |
        has_baseline_creatinine_flow == 0
    ),
    
    no_eligible_intervals = (
      is.na(any_eligible_interval) |
        any_eligible_interval == 0
    ),
    
    no_longitudinal_rows = (
      is.na(n_intervals) |
        n_intervals == 0
    ),
    
    no_final_weights = (
      is.na(has_weight) |
        has_weight == 0
    ),
    
    no_lagged_map = (
      is.na(has_map_lag) |
        has_map_lag == 0
    ),
    
    no_lagged_hr = (
      is.na(has_hr_lag) |
        has_hr_lag == 0
    ),
    
    no_lagged_fluid = (
      is.na(has_fluid_lag) |
        has_fluid_lag == 0
    ),
    
    no_lagged_vaso = (
      is.na(has_vaso_lag) |
        has_vaso_lag == 0
    )
  )

# -----------------------------------------------------------------------------
# Mutually exclusive exclusion reasons
# -----------------------------------------------------------------------------

breakdown <- candidate %>%
  dplyr::mutate(
    exclusion_reason = dplyr::case_when(
      included_long24r ~
        "Included in revised longitudinal analytic dataset",
      
      missing_postop_creatinine ~
        "No evaluable postoperative creatinine within 7 days",
      
      missing_baseline_creatinine ~
        "No baseline creatinine available",
      
      no_eligible_intervals ~
        "No eligible intervals after longitudinal restructuring",
      
      no_longitudinal_rows ~
        "No longitudinal interval rows available",
      
      no_final_weights ~
        "No estimable final stabilised weights",
      
      no_lagged_map ~
        "No lagged MAP information",
      
      no_lagged_hr ~
        "No lagged heart rate information",
      
      no_lagged_fluid ~
        "No lagged fluid exposure information",
      
      no_lagged_vaso ~
        "No lagged vasopressor exposure information",
      
      TRUE ~
        "Other unresolved discrepancy"
    )
  ) %>%
  dplyr::count(exclusion_reason, name = "procedures") %>%
  dplyr::mutate(
    denominator = nrow(candidate),
    percent = 100 * procedures / denominator
  ) %>%
  dplyr::arrange(dplyr::desc(procedures))

# -----------------------------------------------------------------------------
# Non-exclusive detail
# -----------------------------------------------------------------------------

detail <- candidate %>%
  dplyr::summarise(
    denominator = dplyr::n(),
    included_long24r = sum(included_long24r, na.rm = TRUE),
    missing_postop_creatinine = sum(missing_postop_creatinine, na.rm = TRUE),
    missing_baseline_creatinine = sum(missing_baseline_creatinine, na.rm = TRUE),
    no_eligible_intervals = sum(no_eligible_intervals, na.rm = TRUE),
    no_longitudinal_rows = sum(no_longitudinal_rows, na.rm = TRUE),
    no_final_weights = sum(no_final_weights, na.rm = TRUE),
    no_lagged_map = sum(no_lagged_map, na.rm = TRUE),
    no_lagged_hr = sum(no_lagged_hr, na.rm = TRUE),
    no_lagged_fluid = sum(no_lagged_fluid, na.rm = TRUE),
    no_lagged_vaso = sum(no_lagged_vaso, na.rm = TRUE)
  )

# -----------------------------------------------------------------------------
# Consistency checks
# -----------------------------------------------------------------------------

check <- tibble::tibble(
  metric = c(
    "duration_ge_120_candidate",
    "included_long24r",
    "excluded_from_long24r"
  ),
  procedures = c(
    nrow(candidate),
    sum(candidate$included_long24r, na.rm = TRUE),
    sum(!candidate$included_long24r, na.rm = TRUE)
  )
)

# -----------------------------------------------------------------------------
# Print
# -----------------------------------------------------------------------------

message("==== Consistency check ====")
print(tibble::as_tibble(check), n = Inf, width = Inf)

message("==== Mutually exclusive breakdown ====")
print(tibble::as_tibble(breakdown), n = Inf, width = Inf)

message("==== Non-exclusive detail ====")
print(tibble::as_tibble(detail), n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  check,
  file.path(TABLES_DIR, "48R_longitudinal_dataset_loss_check.csv")
)

readr::write_csv(
  breakdown,
  file.path(TABLES_DIR, "48R_longitudinal_dataset_loss_breakdown.csv")
)

readr::write_csv(
  detail,
  file.path(TABLES_DIR, "48R_longitudinal_dataset_loss_detail.csv")
)

saveRDS(
  breakdown,
  file.path(DERIVED_DIR, "48R_longitudinal_dataset_loss_breakdown.rds")
)

saveRDS(
  detail,
  file.path(DERIVED_DIR, "48R_longitudinal_dataset_loss_detail.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "48R_longitudinal_dataset_loss_breakdown.csv"))
message(file.path(TABLES_DIR, "48R_longitudinal_dataset_loss_detail.csv"))
message("==== AUDIT COMPLETE ====")