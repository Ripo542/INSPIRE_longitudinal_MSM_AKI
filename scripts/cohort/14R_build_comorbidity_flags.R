# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 14R_build_comorbidity_flags.R
# Purpose:
#   Build simple baseline comorbidity flags from ICD-10-CM diagnosis codes.
#   Flags are defined using diagnoses recorded before anaesthesia start (t0).
#
#   Intended use:
#     - baseline risk adjustment
#     - not a full Charlson/Elixhauser reconstruction
# =============================================================================

SCRIPT <- "14R_build_comorbidity_flags"

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
  library(janitor)
  library(stringr)
})

message("==== BUILD BASELINE COMORBIDITY FLAGS ====")

# -----------------------------------------------------------------------------
# LOAD BASE COHORT
# -----------------------------------------------------------------------------

base <- readRDS(file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds")) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id),
    t0 = as.numeric(t0)
  ) %>%
  dplyr::filter(has_valid_window)

message("Base valid procedures: ", dplyr::n_distinct(base$op_id))

# -----------------------------------------------------------------------------
# LOAD DIAGNOSIS
# -----------------------------------------------------------------------------

diag_path <- file.path(RAW_DIR, "diagnosis.csv")
stopifnot(file.exists(diag_path))

diag <- data.table::fread(
  diag_path,
  select = c("subject_id", "chart_time", "icd10_cm"),
  showProgress = TRUE
) %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    subject_id = as.character(subject_id),
    chart_time = as.numeric(chart_time),
    icd10_cm = toupper(trimws(as.character(icd10_cm))),
    icd3 = substr(icd10_cm, 1, 3)
  ) %>%
  dplyr::filter(!is.na(icd10_cm), icd10_cm != "")

message("Diagnosis rows loaded: ", nrow(diag))
message("Diagnosis subjects: ", dplyr::n_distinct(diag$subject_id))

# -----------------------------------------------------------------------------
# TEMPORAL LINK TO PROCEDURES
# -----------------------------------------------------------------------------
# Only diagnoses recorded before or at anaesthesia start are considered baseline.
# A broad 5-year lookback is used to avoid remote irrelevant coding while retaining
# chronic disease diagnoses.

lookback_days <- 365 * 5

diag_op <- base %>%
  dplyr::select(op_id, subject_id, t0) %>%
  dplyr::inner_join(diag, by = "subject_id") %>%
  dplyr::mutate(
    rel_to_t0_min = chart_time - t0
  ) %>%
  dplyr::filter(
    rel_to_t0_min <= 0,
    rel_to_t0_min >= -lookback_days * 24 * 60
  )

message("Operation-diagnosis rows within lookback: ", nrow(diag_op))
message("Procedures with any baseline diagnosis: ", dplyr::n_distinct(diag_op$op_id))

# -----------------------------------------------------------------------------
# ICD-10-CM FLAGS
# -----------------------------------------------------------------------------
# Conservative simple flags:
#   diabetes: E10, E11, E13, E14
#   chronic kidney disease: N18
#   hypertension: I10
#   ischaemic heart disease: I20-I25
#   atrial fibrillation/flutter: I48
#   stroke/cerebrovascular disease: I60-I69
#   chronic pulmonary disease: J40-J47
#
# These are deliberately broad baseline flags, not severity-coded comorbidity scores.

diag_flagged <- diag_op %>%
  dplyr::mutate(
    diabetes = as.integer(icd3 %in% c("E10", "E11", "E13", "E14")),
    ckd = as.integer(icd3 == "N18"),
    hypertension = as.integer(icd3 == "I10"),
    ischaemic_heart_disease = as.integer(icd3 >= "I20" & icd3 <= "I25"),
    atrial_fibrillation = as.integer(icd3 == "I48"),
    cerebrovascular_disease = as.integer(icd3 >= "I60" & icd3 <= "I69"),
    chronic_pulmonary_disease = as.integer(icd3 >= "J40" & icd3 <= "J47")
  )

comorbidity_flags <- diag_flagged %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    diabetes = as.integer(any(diabetes == 1, na.rm = TRUE)),
    ckd = as.integer(any(ckd == 1, na.rm = TRUE)),
    hypertension = as.integer(any(hypertension == 1, na.rm = TRUE)),
    ischaemic_heart_disease = as.integer(any(ischaemic_heart_disease == 1, na.rm = TRUE)),
    atrial_fibrillation = as.integer(any(atrial_fibrillation == 1, na.rm = TRUE)),
    cerebrovascular_disease = as.integer(any(cerebrovascular_disease == 1, na.rm = TRUE)),
    chronic_pulmonary_disease = as.integer(any(chronic_pulmonary_disease == 1, na.rm = TRUE)),
    n_baseline_diagnosis_codes = dplyr::n_distinct(icd10_cm),
    .groups = "drop"
  )

# Add zeroes for procedures without baseline diagnosis codes
comorbidity_flags <- base %>%
  dplyr::select(op_id, subject_id) %>%
  dplyr::left_join(comorbidity_flags, by = "op_id") %>%
  dplyr::mutate(
    dplyr::across(
      c(
        diabetes,
        ckd,
        hypertension,
        ischaemic_heart_disease,
        atrial_fibrillation,
        cerebrovascular_disease,
        chronic_pulmonary_disease
      ),
      ~ dplyr::coalesce(.x, 0L)
    ),
    n_baseline_diagnosis_codes = dplyr::coalesce(n_baseline_diagnosis_codes, 0L),
    any_baseline_diagnosis = as.integer(n_baseline_diagnosis_codes > 0)
  )

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

qc_flags <- comorbidity_flags %>%
  dplyr::summarise(
    n_procedures = dplyr::n(),
    n_unique_patients = dplyr::n_distinct(subject_id),
    
    n_any_baseline_diagnosis = sum(any_baseline_diagnosis == 1, na.rm = TRUE),
    any_baseline_diagnosis_pct = 100 * mean(any_baseline_diagnosis == 1, na.rm = TRUE),
    
    n_diabetes = sum(diabetes == 1, na.rm = TRUE),
    diabetes_pct = 100 * mean(diabetes == 1, na.rm = TRUE),
    
    n_ckd = sum(ckd == 1, na.rm = TRUE),
    ckd_pct = 100 * mean(ckd == 1, na.rm = TRUE),
    
    n_hypertension = sum(hypertension == 1, na.rm = TRUE),
    hypertension_pct = 100 * mean(hypertension == 1, na.rm = TRUE),
    
    n_ischaemic_heart_disease = sum(ischaemic_heart_disease == 1, na.rm = TRUE),
    ischaemic_heart_disease_pct = 100 * mean(ischaemic_heart_disease == 1, na.rm = TRUE),
    
    n_atrial_fibrillation = sum(atrial_fibrillation == 1, na.rm = TRUE),
    atrial_fibrillation_pct = 100 * mean(atrial_fibrillation == 1, na.rm = TRUE),
    
    n_cerebrovascular_disease = sum(cerebrovascular_disease == 1, na.rm = TRUE),
    cerebrovascular_disease_pct = 100 * mean(cerebrovascular_disease == 1, na.rm = TRUE),
    
    n_chronic_pulmonary_disease = sum(chronic_pulmonary_disease == 1, na.rm = TRUE),
    chronic_pulmonary_disease_pct = 100 * mean(chronic_pulmonary_disease == 1, na.rm = TRUE),
    
    diagnosis_codes_median = stats::median(n_baseline_diagnosis_codes, na.rm = TRUE),
    diagnosis_codes_p25 = as.numeric(stats::quantile(n_baseline_diagnosis_codes, 0.25, na.rm = TRUE)),
    diagnosis_codes_p75 = as.numeric(stats::quantile(n_baseline_diagnosis_codes, 0.75, na.rm = TRUE))
  )

message("==== Comorbidity flag QC ====")
print(tibble::as_tibble(qc_flags), width = Inf)

# Top baseline codes after temporal restriction
top_codes <- diag_op %>%
  dplyr::count(icd10_cm, sort = TRUE) %>%
  dplyr::slice_head(n = 50)

message("==== Top baseline diagnosis codes within lookback ====")
print(tibble::as_tibble(top_codes), n = 50)

# QC by flag: top codes contributing
flag_code_qc <- diag_flagged %>%
  dplyr::mutate(
    flag = dplyr::case_when(
      diabetes == 1 ~ "diabetes",
      ckd == 1 ~ "ckd",
      hypertension == 1 ~ "hypertension",
      ischaemic_heart_disease == 1 ~ "ischaemic_heart_disease",
      atrial_fibrillation == 1 ~ "atrial_fibrillation",
      cerebrovascular_disease == 1 ~ "cerebrovascular_disease",
      chronic_pulmonary_disease == 1 ~ "chronic_pulmonary_disease",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(flag)) %>%
  dplyr::count(flag, icd10_cm, sort = TRUE) %>%
  dplyr::group_by(flag) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::ungroup()

message("==== Top contributing codes per comorbidity flag ====")
print(tibble::as_tibble(flag_code_qc), n = 100)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  comorbidity_flags,
  file.path(DERIVED_DIR, "14R_baseline_comorbidity_flags.rds")
)

readr::write_csv(
  comorbidity_flags,
  file.path(TABLES_DIR, "14R_baseline_comorbidity_flags.csv")
)

readr::write_csv(
  qc_flags,
  file.path(QC_DIR, "14R_comorbidity_flags_qc.csv")
)

readr::write_csv(
  top_codes,
  file.path(QC_DIR, "14R_top_baseline_diagnosis_codes.csv")
)

readr::write_csv(
  flag_code_qc,
  file.path(QC_DIR, "14R_flag_code_qc.csv")
)

write_qc_json(
  list(
    script = SCRIPT,
    timestamp = as.character(Sys.time()),
    lookback_days = lookback_days,
    comorbidity_qc = qc_flags,
    top_baseline_codes = top_codes,
    flag_code_qc = flag_code_qc
  ),
  "14R_build_comorbidity_flags.json"
)

message("==== COMORBIDITY FLAGS COMPLETE ====")
message("Saved: ", file.path(DERIVED_DIR, "14R_baseline_comorbidity_flags.rds"))