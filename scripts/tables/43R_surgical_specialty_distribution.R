# =============================================================================
# 43R_surgical_specialty_distribution.R
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

abd <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_elective_abdominal_major.rds"))
full <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_full_major_noncardiac.rds"))

fmt <- function(n, total) {
  sprintf("%d (%.1f%%)", n, 100 * n / total)
}

# Abdominal
abd_counts <- abd %>%
  dplyr::count(department_final) %>%
  dplyr::mutate(val = fmt(n, sum(n))) %>%
  dplyr::select(department_final, abd = val)

# Full
full_counts <- full %>%
  dplyr::count(department_final) %>%
  dplyr::mutate(val = fmt(n, sum(n))) %>%
  dplyr::select(department_final, full = val)

# Merge
tab <- full_counts %>%
  dplyr::left_join(abd_counts, by = "department_final") %>%
  dplyr::mutate(
    abd = ifelse(is.na(abd), "–", abd)
  ) %>%
  dplyr::rename(Surgical_specialty = department_final)

message("==== Surgical specialty distribution ====")
print(tab, n = Inf)

# Save
readr::write_csv(
  tab,
  file.path(TABLES_DIR, "43R_surgical_specialty_distribution.csv")
)

message("Saved:")
message(file.path(TABLES_DIR, "43R_surgical_specialty_distribution.csv"))