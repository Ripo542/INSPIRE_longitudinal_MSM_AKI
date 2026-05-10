# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 41R_table1_baseline_characteristics.R
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== BUILD TABLE 1: REVISED COHORT CHARACTERISTICS ====")

abd <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_elective_abdominal_major.rds"))
full <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_full_major_noncardiac.rds"))

fmt_n_pct <- function(x, denom, digits = 1) {
  sprintf(paste0("%d (%.", digits, "f%%)"), sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
}

fmt_mean_sd <- function(x, digits = 1) {
  sprintf(
    paste0("%.", digits, "f ± %.", digits, "f"),
    mean(x, na.rm = TRUE),
    stats::sd(x, na.rm = TRUE)
  )
}

fmt_median_iqr <- function(x, digits = 1) {
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f–%.", digits, "f]"),
    q[[2]], q[[1]], q[[3]]
  )
}

summarise_table1 <- function(df) {
  
  n_proc <- nrow(df)
  
  tibble::tibble(
    Characteristic = c(
      "Surgical procedures, n",
      "Unique patients, n",
      "Age, yr",
      "Female sex",
      "Body mass index, kg m−2",
      "Baseline creatinine, mg dl−1",
      "Baseline haemoglobin, g dl−1",
      "Baseline albumin, g dl−1",
      "ASA physical status I",
      "ASA physical status II",
      "ASA physical status III",
      "ASA physical status IV–V",
      "Diabetes",
      "Chronic kidney disease",
      "Ischaemic heart disease",
      "Chronic pulmonary disease",
      "Procedure duration, min",
      "5-min intervals per procedure",
      "Any MAP <65 mmHg",
      "Hypotension duration, min",
      "Any fluid administration",
      "Fluid exposure, 5-min intervals",
      "Any vasopressor administration",
      "Vasopressor exposure, 5-min intervals",
      "Postoperative AKI"
    ),
    
    Value = c(
      as.character(n_proc),
      as.character(dplyr::n_distinct(df$subject_id)),
      fmt_mean_sd(df$age, 1),
      fmt_n_pct(df$female == 1, n_proc, 1),
      fmt_median_iqr(df$bmi, 1),
      fmt_median_iqr(df$baseline_creatinine, 2),
      fmt_median_iqr(df$baseline_haemoglobin, 1),
      fmt_median_iqr(df$baseline_albumin, 1),
      
      fmt_n_pct(df$asa_final == "I", n_proc, 1),
      fmt_n_pct(df$asa_final == "II", n_proc, 1),
      fmt_n_pct(df$asa_final == "III", n_proc, 1),
      fmt_n_pct(df$asa_final == "IV_V", n_proc, 1),
      
      fmt_n_pct(df$diabetes == 1, n_proc, 1),
      fmt_n_pct(df$ckd == 1, n_proc, 1),
      fmt_n_pct(df$ischaemic_heart_disease == 1, n_proc, 1),
      fmt_n_pct(df$chronic_pulmonary_disease == 1, n_proc, 1),
      
      fmt_median_iqr(df$n_intervals * 5, 0),
      fmt_median_iqr(df$n_intervals, 0),
      
      fmt_n_pct(df$hypotension_time_10min > 0, n_proc, 1),
      fmt_median_iqr(df$hypotension_time_10min * 10, 0),
      
      fmt_n_pct(df$fluid_intervals_per_5 > 0, n_proc, 1),
      fmt_median_iqr(df$fluid_intervals_per_5, 0),
      
      fmt_n_pct(df$vasopressor_intervals_per_5 > 0, n_proc, 1),
      fmt_median_iqr(df$vasopressor_intervals_per_5, 0),
      
      fmt_n_pct(df$aki_creat_7d == 1, n_proc, 1)
    )
  )
}

tab_abd <- summarise_table1(abd) %>%
  dplyr::rename(`Primary abdominal cohort` = Value)

tab_full <- summarise_table1(full) %>%
  dplyr::rename(`Full major non-cardiac cohort` = Value)

table1 <- tab_abd %>%
  dplyr::left_join(tab_full, by = "Characteristic")

message("==== Table 1 ====")
print(table1, n = Inf)

readr::write_csv(
  table1,
  file.path(TABLES_DIR, "41R_table1_revised_cohort_characteristics.csv")
)

saveRDS(
  table1,
  file.path(DERIVED_DIR, "41R_table1_revised_cohort_characteristics.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "41R_table1_revised_cohort_characteristics.csv"))
message("==== TABLE 1 COMPLETE ====")