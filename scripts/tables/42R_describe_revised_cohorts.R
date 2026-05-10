# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 42R_describe_revised_cohorts.R
# Purpose:
#   Describe revised primary abdominal and full sensitivity cohorts:
#   - procedures / patients
#   - intervals
#   - MAP / HR missingness
#   - hypotension exposure
#   - fluid / vasopressor exposure
#   - AKI incidence
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== DESCRIBE REVISED COHORTS ====")

# -----------------------------------------------------------------------------
# Load revised longitudinal dataset
# -----------------------------------------------------------------------------

long <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds"))

message("Longitudinal rows: ", nrow(long))
message("Operations: ", dplyr::n_distinct(long$op_id))

# -----------------------------------------------------------------------------
# Helper
# -----------------------------------------------------------------------------

summarise_cohort <- function(data, flag, label) {
  
  d <- data %>%
    dplyr::filter(.data[[flag]] == 1)
  
  op_level <- d %>%
    dplyr::group_by(op_id) %>%
    dplyr::summarise(
      subject_id = dplyr::first(subject_id),
      aki_creat_7d = dplyr::first(aki_creat_7d),
      n_intervals = dplyr::n(),
      any_map = any(!is.na(map_lag1)),
      any_hr = any(!is.na(hr_lag1)),
      any_hypotension = any(hypotension65_lag1 == 1, na.rm = TRUE),
      any_fluid = any(fluid_any == 1, na.rm = TRUE),
      any_vasopressor = any(vasopressor_any == 1, na.rm = TRUE),
      .groups = "drop"
    )
  
  tibble::tibble(
    cohort = label,
    
    n_procedures = dplyr::n_distinct(d$op_id),
    n_patients = dplyr::n_distinct(d$subject_id),
    n_intervals = nrow(d),
    
    intervals_median = stats::median(op_level$n_intervals, na.rm = TRUE),
    intervals_p25 = as.numeric(stats::quantile(op_level$n_intervals, 0.25, na.rm = TRUE)),
    intervals_p75 = as.numeric(stats::quantile(op_level$n_intervals, 0.75, na.rm = TRUE)),
    
    duration_median_min = stats::median(op_level$n_intervals * 5, na.rm = TRUE),
    duration_p25_min = as.numeric(stats::quantile(op_level$n_intervals * 5, 0.25, na.rm = TRUE)),
    duration_p75_min = as.numeric(stats::quantile(op_level$n_intervals * 5, 0.75, na.rm = TRUE)),
    
    map_missing_pct = 100 * mean(is.na(d$map_lag1)),
    hr_missing_pct = 100 * mean(is.na(d$hr_lag1)),
    
    hypotension_interval_pct = 100 * mean(d$hypotension65_lag1 == 1, na.rm = TRUE),
    fluid_interval_pct = 100 * mean(d$fluid_any == 1, na.rm = TRUE),
    vasopressor_interval_pct = 100 * mean(d$vasopressor_any == 1, na.rm = TRUE),
    
    procedures_any_hypotension = sum(op_level$any_hypotension, na.rm = TRUE),
    procedures_any_hypotension_pct = 100 * mean(op_level$any_hypotension, na.rm = TRUE),
    
    procedures_any_fluid = sum(op_level$any_fluid, na.rm = TRUE),
    procedures_any_fluid_pct = 100 * mean(op_level$any_fluid, na.rm = TRUE),
    
    procedures_any_vasopressor = sum(op_level$any_vasopressor, na.rm = TRUE),
    procedures_any_vasopressor_pct = 100 * mean(op_level$any_vasopressor, na.rm = TRUE),
    
    aki_events = sum(op_level$aki_creat_7d == 1, na.rm = TRUE),
    aki_risk_pct = 100 * mean(op_level$aki_creat_7d == 1, na.rm = TRUE)
  )
}

# -----------------------------------------------------------------------------
# Cohort summaries
# -----------------------------------------------------------------------------

cohort_summary <- dplyr::bind_rows(
  summarise_cohort(
    long,
    "revised_primary_abdominal_complete",
    "Primary: elective major abdominal surgery"
  ),
  summarise_cohort(
    long,
    "revised_primary_full_complete",
    "Sensitivity: full major non-cardiac surgery"
  )
)

message("==== Revised cohort summary ====")
print(cohort_summary, width = Inf)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  cohort_summary,
  file.path(TABLES_DIR, "42R_revised_cohort_summary.csv")
)

saveRDS(
  cohort_summary,
  file.path(DERIVED_DIR, "42R_revised_cohort_summary.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "42R_revised_cohort_summary.csv"))
message("==== REVISED COHORT DESCRIPTION COMPLETE ====")