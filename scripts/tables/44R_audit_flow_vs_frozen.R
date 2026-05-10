# =============================================================================
# 44R_audit_flow_vs_frozen.R
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

flow <- readRDS(file.path(DERIVED_DIR, "44R_flowchart_counts.rds"))

long <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds"))

full_frozen <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_full_major_noncardiac.rds")) %>%
  mutate(op_id = as.character(op_id))

abd_frozen <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_elective_abdominal_major.rds")) %>%
  mutate(op_id = as.character(op_id))

# Recreate upstream flags from revised longitudinal object
op_long <- long %>%
  group_by(op_id) %>%
  summarise(
    subject_id = first(subject_id),
    duration_min = first(duration_min),
    revised_primary_full_complete = first(revised_primary_full_complete),
    revised_primary_abdominal_complete = first(revised_primary_abdominal_complete),
    any_eligible_interval = n() > 0,
    has_weight = all(!is.na(final_sw_joint_trunc_1_99)),
    n_intervals = n(),
    .groups = "drop"
  )

audit <- tibble(
  comparison = c(
    "full_long_flag",
    "full_frozen",
    "abd_long_flag",
    "abd_frozen"
  ),
  n_ops = c(
    n_distinct(op_long$op_id[op_long$revised_primary_full_complete == 1]),
    n_distinct(full_frozen$op_id),
    n_distinct(op_long$op_id[op_long$revised_primary_abdominal_complete == 1]),
    n_distinct(abd_frozen$op_id)
  )
)

print(audit)

full_missing <- setdiff(
  op_long$op_id[op_long$revised_primary_full_complete == 1],
  full_frozen$op_id
)

abd_missing <- setdiff(
  op_long$op_id[op_long$revised_primary_abdominal_complete == 1],
  abd_frozen$op_id
)

cat("\nFull missing from frozen:", length(full_missing), "\n")
cat("Abdominal missing from frozen:", length(abd_missing), "\n")

missing_summary <- op_long %>%
  filter(op_id %in% c(full_missing, abd_missing)) %>%
  summarise(
    n_ops = n(),
    n_intervals_median = median(n_intervals, na.rm = TRUE),
    n_intervals_min = min(n_intervals, na.rm = TRUE),
    n_intervals_p25 = quantile(n_intervals, 0.25, na.rm = TRUE),
    n_intervals_p75 = quantile(n_intervals, 0.75, na.rm = TRUE),
    missing_weight = sum(!has_weight, na.rm = TRUE)
  )

print(missing_summary)

write_csv(audit, file.path(QC_DIR, "44R_audit_flow_vs_frozen_counts.csv"))
write_csv(missing_summary, file.path(QC_DIR, "44R_audit_flow_vs_frozen_missing_summary.csv"))