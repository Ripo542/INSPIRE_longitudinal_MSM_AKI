# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 35R_weight_diagnostics.R
# Purpose:
#   Diagnose final stabilised weights for revised MSM cohorts.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== WEIGHT DIAGNOSTICS: REVISED COHORTS ====")

long <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds"))

summarise_weights <- function(data, flag, label) {
  
  w <- data %>%
    dplyr::filter(.data[[flag]] == 1) %>%
    dplyr::arrange(op_id, interval_id) %>%
    dplyr::group_by(op_id) %>%
    dplyr::summarise(
      final_sw_joint = dplyr::last(final_sw_joint),
      final_sw_joint_trunc_1_99 = dplyr::last(final_sw_joint_trunc_1_99),
      final_sw_joint_trunc_0_5_99_5 = dplyr::last(final_sw_joint_trunc_0_5_99_5),
      .groups = "drop"
    )
  
  tibble::tibble(
    cohort = label,
    n_ops = nrow(w),
    
    untruncated_mean = mean(w$final_sw_joint, na.rm = TRUE),
    untruncated_median = stats::median(w$final_sw_joint, na.rm = TRUE),
    untruncated_p95 = as.numeric(stats::quantile(w$final_sw_joint, 0.95, na.rm = TRUE)),
    untruncated_p99 = as.numeric(stats::quantile(w$final_sw_joint, 0.99, na.rm = TRUE)),
    untruncated_max = max(w$final_sw_joint, na.rm = TRUE),
    
    trunc_1_99_mean = mean(w$final_sw_joint_trunc_1_99, na.rm = TRUE),
    trunc_1_99_median = stats::median(w$final_sw_joint_trunc_1_99, na.rm = TRUE),
    trunc_1_99_p95 = as.numeric(stats::quantile(w$final_sw_joint_trunc_1_99, 0.95, na.rm = TRUE)),
    trunc_1_99_p99 = as.numeric(stats::quantile(w$final_sw_joint_trunc_1_99, 0.99, na.rm = TRUE)),
    trunc_1_99_max = max(w$final_sw_joint_trunc_1_99, na.rm = TRUE),
    
    trunc_0_5_99_5_mean = mean(w$final_sw_joint_trunc_0_5_99_5, na.rm = TRUE),
    trunc_0_5_99_5_median = stats::median(w$final_sw_joint_trunc_0_5_99_5, na.rm = TRUE),
    trunc_0_5_99_5_p95 = as.numeric(stats::quantile(w$final_sw_joint_trunc_0_5_99_5, 0.95, na.rm = TRUE)),
    trunc_0_5_99_5_p995 = as.numeric(stats::quantile(w$final_sw_joint_trunc_0_5_99_5, 0.995, na.rm = TRUE)),
    trunc_0_5_99_5_max = max(w$final_sw_joint_trunc_0_5_99_5, na.rm = TRUE)
  )
}

weight_diag <- dplyr::bind_rows(
  summarise_weights(
    long,
    "revised_primary_abdominal_complete",
    "Primary abdominal cohort"
  ),
  summarise_weights(
    long,
    "revised_primary_full_complete",
    "Full major non-cardiac cohort"
  )
)

message("==== Weight diagnostics ====")
print(weight_diag, width = Inf)

readr::write_csv(
  weight_diag,
  file.path(TABLES_DIR, "35R_weight_diagnostics_revised_cohorts.csv")
)

saveRDS(
  weight_diag,
  file.path(DERIVED_DIR, "35R_weight_diagnostics_revised_cohorts.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "35R_weight_diagnostics_revised_cohorts.csv"))
message("==== WEIGHT DIAGNOSTICS COMPLETE ====")