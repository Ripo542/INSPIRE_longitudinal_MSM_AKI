# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 49R_table_missing_included_vs_excluded.R
# Purpose:
#   Compare included vs excluded procedures after duration >=120 min and report
#   missingness by variable.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== INCLUDED VS EXCLUDED TABLE WITH MISSINGNESS ====")

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
    baseline_creatinine,
    has_postop_creatinine_7d,
    aki_creat_7d
  )

labs <- readRDS(file.path(DERIVED_DIR, "12R_baseline_haemoglobin_albumin.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id)) %>%
  dplyr::select(
    op_id,
    baseline_haemoglobin,
    baseline_albumin
  )

long_final <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

included_ids <- unique(long_final$op_id)

# -----------------------------------------------------------------------------
# Candidate cohort: duration >=120 min
# -----------------------------------------------------------------------------

df <- base %>%
  dplyr::left_join(aki, by = "op_id") %>%
  dplyr::left_join(labs, by = "op_id") %>%
  dplyr::mutate(
    duration_ge_120 = intraop_duration_min >= 120,
    included = op_id %in% included_ids,
    
    female_final = dplyr::case_when(
      "female" %in% names(.) ~ as.integer(female),
      TRUE ~ as.integer(sex_raw == "F")
    ),
    
    department_final = dplyr::coalesce(
      as.character(department_clean),
      as.character(department)
    )
  ) %>%
  dplyr::filter(duration_ge_120 == TRUE)

message("Candidate procedures >=120 min: ", nrow(df))
message("Included in revised longitudinal dataset: ", sum(df$included))
message("Excluded before revised longitudinal dataset: ", sum(!df$included))

# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------

fmt_n_pct <- function(x) {
  sprintf(
    "%d (%.1f%%)",
    sum(x, na.rm = TRUE),
    100 * mean(x, na.rm = TRUE)
  )
}

fmt_mean_sd <- function(x, digits = 1) {
  sprintf(
    paste0("%.", digits, "f ± %.", digits, "f"),
    mean(x, na.rm = TRUE),
    stats::sd(x, na.rm = TRUE)
  )
}

fmt_median_iqr <- function(x, digits = 1) {
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f–%.", digits, "f]"),
    q[[2]], q[[1]], q[[3]]
  )
}

fmt_missing <- function(x) {
  sprintf(
    "%d (%.1f%%)",
    sum(is.na(x)),
    100 * mean(is.na(x))
  )
}

# -----------------------------------------------------------------------------
# Summary table by inclusion status
# -----------------------------------------------------------------------------

summarise_group <- function(data) {
  tibble::tibble(
    `Procedures, n` = as.character(nrow(data)),
    `Unique patients, n` = as.character(dplyr::n_distinct(data$subject_id)),
    
    `Age, yr` = fmt_mean_sd(data$age, 1),
    `Female sex` = fmt_n_pct(data$female_final == 1),
    `Body mass index, kg m−2` = fmt_median_iqr(data$bmi, 1),
    
    `ASA physical status I` = fmt_n_pct(as.character(data$asa) %in% c("1", "I")),
    `ASA physical status II` = fmt_n_pct(as.character(data$asa) %in% c("2", "II")),
    `ASA physical status III` = fmt_n_pct(as.character(data$asa) %in% c("3", "III")),
    `ASA physical status IV–V` = fmt_n_pct(as.character(data$asa) %in% c("4", "5", "IV", "V")),
    
    `Baseline creatinine, mg dl−1` = fmt_median_iqr(data$baseline_creatinine, 2),
    `Baseline haemoglobin, g dl−1` = fmt_median_iqr(data$baseline_haemoglobin, 1),
    `Baseline albumin, g dl−1` = fmt_median_iqr(data$baseline_albumin, 1),
    
    `Intraoperative duration, min` = fmt_median_iqr(data$intraop_duration_min, 0),
    `Postoperative AKI evaluable` = fmt_n_pct(!is.na(data$aki_creat_7d)),
    `Postoperative AKI` = fmt_n_pct(data$aki_creat_7d == 1)
  )
}

included_tbl <- summarise_group(df %>% dplyr::filter(included == TRUE)) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Characteristic",
    values_to = "Included"
  )

excluded_tbl <- summarise_group(df %>% dplyr::filter(included == FALSE)) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Characteristic",
    values_to = "Excluded"
  )

table_included_excluded <- included_tbl %>%
  dplyr::left_join(excluded_tbl, by = "Characteristic")

# -----------------------------------------------------------------------------
# Missingness table by inclusion status
# -----------------------------------------------------------------------------

missingness_group <- function(data) {
  tibble::tibble(
    `Age` = fmt_missing(data$age),
    `Sex` = fmt_missing(data$female_final),
    `ASA physical status` = fmt_missing(data$asa),
    `Body mass index` = fmt_missing(data$bmi),
    `Baseline creatinine` = fmt_missing(data$baseline_creatinine),
    `Baseline haemoglobin` = fmt_missing(data$baseline_haemoglobin),
    `Baseline albumin` = fmt_missing(data$baseline_albumin),
    `Surgical specialty` = fmt_missing(data$department_final),
    `Postoperative AKI outcome` = fmt_missing(data$aki_creat_7d)
  )
}

missing_included <- missingness_group(df %>% dplyr::filter(included == TRUE)) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Included_missing"
  )

missing_excluded <- missingness_group(df %>% dplyr::filter(included == FALSE)) %>%
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "Variable",
    values_to = "Excluded_missing"
  )

table_missingness <- missing_included %>%
  dplyr::left_join(missing_excluded, by = "Variable")

# -----------------------------------------------------------------------------
# Print
# -----------------------------------------------------------------------------

message("==== Included vs excluded characteristics ====")
print(tibble::as_tibble(table_included_excluded), n = Inf, width = Inf)

message("==== Missingness by variable ====")
print(tibble::as_tibble(table_missingness), n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  table_included_excluded,
  file.path(TABLES_DIR, "49R_included_vs_excluded_characteristics.csv")
)

readr::write_csv(
  table_missingness,
  file.path(TABLES_DIR, "49R_included_vs_excluded_missingness.csv")
)

saveRDS(
  table_included_excluded,
  file.path(DERIVED_DIR, "49R_included_vs_excluded_characteristics.rds")
)

saveRDS(
  table_missingness,
  file.path(DERIVED_DIR, "49R_included_vs_excluded_missingness.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "49R_included_vs_excluded_characteristics.csv"))
message(file.path(TABLES_DIR, "49R_included_vs_excluded_missingness.csv"))
message("==== INCLUDED VS EXCLUDED COMPLETE ====")