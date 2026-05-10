# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 11N_build_aki_outcome.R
# Purpose:
#   Build creatinine-based postoperative AKI outcome within 7 days after surgery.
#   Uses INSPIRE labs.csv + operation time windows.
# =============================================================================

SCRIPT <- "11N_build_aki_outcome"

source("R/00_setup/00N_setup_project.R")

message("==== BUILD AKI OUTCOME ====")

# -----------------------------------------------------------------------------
# LOAD BASE COHORT
# -----------------------------------------------------------------------------

base <- readRDS(file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds")) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id),
    t0 = as.numeric(t0),
    t1 = as.numeric(t1)
  ) %>%
  dplyr::filter(has_valid_window)

message("Base valid operations: ", dplyr::n_distinct(base$op_id))

# -----------------------------------------------------------------------------
# LOAD LABS
# -----------------------------------------------------------------------------

labs_path <- file.path(RAW_DIR, "labs.csv")
stopifnot(file.exists(labs_path))

labs <- data.table::fread(
  labs_path,
  select = c("subject_id", "chart_time", "item_name", "value"),
  showProgress = TRUE
) %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    subject_id = as.character(subject_id),
    chart_time = as.numeric(chart_time),
    value = as.numeric(value),
    item_name = as.character(item_name)
  )

message("Labs rows loaded: ", nrow(labs))
message("Lab item names:")
print(sort(unique(labs$item_name)))

# -----------------------------------------------------------------------------
# CREATININE EXTRACTION
# -----------------------------------------------------------------------------

creat_items <- "creatinine"

message("Candidate creatinine items: ", paste(creat_items, collapse = ", "))

if (!"creatinine" %in% unique(labs$item_name)) {
  stop("Creatinine item not found in labs.csv.")
}

creat <- labs %>%
  dplyr::filter(item_name %in% creat_items) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::filter(value > 0)

message("Creatinine rows: ", nrow(creat))
message("Subjects with creatinine: ", dplyr::n_distinct(creat$subject_id))

creat_summary <- creat %>%
  dplyr::group_by(item_name) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_subjects = dplyr::n_distinct(subject_id),
    value_min = min(value, na.rm = TRUE),
    value_p01 = as.numeric(stats::quantile(value, 0.01, na.rm = TRUE)),
    value_median = stats::median(value, na.rm = TRUE),
    value_p99 = as.numeric(stats::quantile(value, 0.99, na.rm = TRUE)),
    value_max = max(value, na.rm = TRUE),
    .groups = "drop"
  )

message("==== Creatinine item summary ====")
print(creat_summary)

readr::write_csv(
  creat_summary,
  file.path(QC_DIR, "11N_creatinine_item_summary.csv")
)

# -----------------------------------------------------------------------------
# JOIN CREATININE TO OPERATIONS BY SUBJECT + TIME
# -----------------------------------------------------------------------------

# Baseline creatinine:
#   primary definition: closest preoperative creatinine within 30 days before t0.
#   fallback: earliest creatinine from t0 to 6h after t0, flagged separately.
#
# Postoperative AKI:
#   maximum creatinine within 7 days after end of anesthesia/surgery window t1.
#   AKI if postoperative creatinine >= 1.5 x baseline OR increase >=0.3 mg/dL.
#
# Unit assumption:
#   INSPIRE creatinine values appear compatible with mg/dL; we audit ranges.

creat_op <- base %>%
  dplyr::select(op_id, subject_id, t0, t1) %>%
  dplyr::inner_join(creat, by = "subject_id") %>%
  dplyr::mutate(
    rel_to_t0_min = chart_time - t0,
    rel_to_t1_min = chart_time - t1
  )

message("Operation-lab creatinine rows after subject join: ", nrow(creat_op))
message("Operations with any creatinine row: ", dplyr::n_distinct(creat_op$op_id))

# -----------------------------------------------------------------------------
# BASELINE CREATININE
# -----------------------------------------------------------------------------

baseline_preop <- creat_op %>%
  dplyr::filter(rel_to_t0_min < 0, rel_to_t0_min >= -30 * 24 * 60) %>%
  dplyr::mutate(abs_time_to_t0 = abs(rel_to_t0_min)) %>%
  dplyr::arrange(op_id, abs_time_to_t0) %>%
  dplyr::group_by(op_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    op_id,
    baseline_creatinine = value,
    baseline_creatinine_time = chart_time,
    baseline_creatinine_source = "preoperative_30d",
    baseline_creatinine_rel_to_t0_min = rel_to_t0_min
  )

baseline_early <- creat_op %>%
  dplyr::filter(rel_to_t0_min >= 0, rel_to_t0_min <= 6 * 60) %>%
  dplyr::arrange(op_id, rel_to_t0_min) %>%
  dplyr::group_by(op_id) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    op_id,
    baseline_creatinine = value,
    baseline_creatinine_time = chart_time,
    baseline_creatinine_source = "early_perioperative_0_6h",
    baseline_creatinine_rel_to_t0_min = rel_to_t0_min
  )

baseline <- base %>%
  dplyr::select(op_id) %>%
  dplyr::left_join(baseline_preop, by = "op_id") %>%
  dplyr::left_join(
    baseline_early,
    by = "op_id",
    suffix = c("", "_early")
  ) %>%
  dplyr::mutate(
    baseline_creatinine = dplyr::coalesce(baseline_creatinine, baseline_creatinine_early),
    baseline_creatinine_time = dplyr::coalesce(baseline_creatinine_time, baseline_creatinine_time_early),
    baseline_creatinine_source = dplyr::coalesce(baseline_creatinine_source, baseline_creatinine_source_early),
    baseline_creatinine_rel_to_t0_min = dplyr::coalesce(
      baseline_creatinine_rel_to_t0_min,
      baseline_creatinine_rel_to_t0_min_early
    )
  ) %>%
  dplyr::select(
    op_id,
    baseline_creatinine,
    baseline_creatinine_time,
    baseline_creatinine_source,
    baseline_creatinine_rel_to_t0_min
  )

# -----------------------------------------------------------------------------
# POSTOPERATIVE CREATININE WITHIN 7 DAYS
# -----------------------------------------------------------------------------

postop <- creat_op %>%
  dplyr::filter(rel_to_t1_min > 0, rel_to_t1_min <= 7 * 24 * 60) %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    postop_creatinine_max_7d = max(value, na.rm = TRUE),
    postop_creatinine_first_7d = value[order(rel_to_t1_min)][1],
    n_postop_creatinine_7d = dplyr::n(),
    first_postop_creatinine_rel_to_t1_min = min(rel_to_t1_min, na.rm = TRUE),
    last_postop_creatinine_rel_to_t1_min = max(rel_to_t1_min, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# OUTCOME DATASET
# -----------------------------------------------------------------------------

aki <- base %>%
  dplyr::select(
    op_id, subject_id, age, female, asa, bmi, department_clean,
    emop, elective, abdominal_department, duration_ge_120,
    t0, t1, intraop_duration_min
  ) %>%
  dplyr::left_join(baseline, by = "op_id") %>%
  dplyr::left_join(postop, by = "op_id") %>%
  dplyr::mutate(
    has_baseline_creatinine = !is.na(baseline_creatinine),
    has_postop_creatinine_7d = !is.na(postop_creatinine_max_7d),
    
    creatinine_abs_increase_7d = postop_creatinine_max_7d - baseline_creatinine,
    creatinine_ratio_7d = postop_creatinine_max_7d / baseline_creatinine,
    
    aki_creat_7d = dplyr::case_when(
      !has_baseline_creatinine | !has_postop_creatinine_7d ~ NA_integer_,
      creatinine_abs_increase_7d >= 0.3 ~ 1L,
      creatinine_ratio_7d >= 1.5 ~ 1L,
      TRUE ~ 0L
    ),
    
    aki_stage_creat_7d = dplyr::case_when(
      is.na(aki_creat_7d) ~ NA_integer_,
      creatinine_ratio_7d >= 3 | postop_creatinine_max_7d >= 4 ~ 3L,
      creatinine_ratio_7d >= 2 ~ 2L,
      aki_creat_7d == 1L ~ 1L,
      TRUE ~ 0L
    )
  )

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

outcome_qc <- aki %>%
  dplyr::summarise(
    n_ops = dplyr::n(),
    baseline_creatinine_available = sum(has_baseline_creatinine, na.rm = TRUE),
    baseline_creatinine_available_pct = 100 * mean(has_baseline_creatinine, na.rm = TRUE),
    baseline_preop_30d = sum(baseline_creatinine_source == "preoperative_30d", na.rm = TRUE),
    baseline_early_0_6h = sum(baseline_creatinine_source == "early_perioperative_0_6h", na.rm = TRUE),
    postop_creatinine_7d_available = sum(has_postop_creatinine_7d, na.rm = TRUE),
    postop_creatinine_7d_available_pct = 100 * mean(has_postop_creatinine_7d, na.rm = TRUE),
    aki_evaluable = sum(!is.na(aki_creat_7d)),
    aki_events = sum(aki_creat_7d == 1, na.rm = TRUE),
    aki_risk_pct = 100 * mean(aki_creat_7d == 1, na.rm = TRUE),
    baseline_creatinine_median = stats::median(baseline_creatinine, na.rm = TRUE),
    postop_creatinine_max_7d_median = stats::median(postop_creatinine_max_7d, na.rm = TRUE)
  )

message("==== AKI outcome QC ====")
print(outcome_qc)

outcome_qc_primary <- aki %>%
  dplyr::filter(duration_ge_120) %>%
  dplyr::summarise(
    n_ops_duration_ge_120 = dplyr::n(),
    aki_evaluable = sum(!is.na(aki_creat_7d)),
    aki_events = sum(aki_creat_7d == 1, na.rm = TRUE),
    aki_risk_pct = 100 * mean(aki_creat_7d == 1, na.rm = TRUE),
    baseline_available_pct = 100 * mean(has_baseline_creatinine, na.rm = TRUE),
    postop_available_pct = 100 * mean(has_postop_creatinine_7d, na.rm = TRUE)
  )

message("==== AKI outcome QC: duration >=120 min ====")
print(outcome_qc_primary)

outcome_qc_abdominal <- aki %>%
  dplyr::filter(
    duration_ge_120,
    abdominal_department == TRUE,
    elective == TRUE | is.na(elective)
  ) %>%
  dplyr::summarise(
    n_ops_abdominal_elective_duration_ge_120 = dplyr::n(),
    aki_evaluable = sum(!is.na(aki_creat_7d)),
    aki_events = sum(aki_creat_7d == 1, na.rm = TRUE),
    aki_risk_pct = 100 * mean(aki_creat_7d == 1, na.rm = TRUE),
    baseline_available_pct = 100 * mean(has_baseline_creatinine, na.rm = TRUE),
    postop_available_pct = 100 * mean(has_postop_creatinine_7d, na.rm = TRUE)
  )

message("==== AKI outcome QC: abdominal/elective/duration >=120 min ====")
print(outcome_qc_abdominal)

baseline_source_qc <- aki %>%
  dplyr::count(baseline_creatinine_source, sort = TRUE)

message("==== Baseline creatinine source QC ====")
print(baseline_source_qc)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  aki,
  file.path(DERIVED_DIR, "11N_aki_outcome_creatinine_7d.rds")
)

readr::write_csv(
  outcome_qc,
  file.path(QC_DIR, "11N_aki_outcome_qc_all.csv")
)

readr::write_csv(
  outcome_qc_primary,
  file.path(QC_DIR, "11N_aki_outcome_qc_duration_ge_120.csv")
)

readr::write_csv(
  outcome_qc_abdominal,
  file.path(QC_DIR, "11N_aki_outcome_qc_abdominal_elective_duration_ge_120.csv")
)

readr::write_csv(
  baseline_source_qc,
  file.path(QC_DIR, "11N_baseline_creatinine_source_qc.csv")
)

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  raw_dir = RAW_DIR,
  creatinine_items = creat_items,
  outcome_qc_all = outcome_qc,
  outcome_qc_duration_ge_120 = outcome_qc_primary,
  outcome_qc_abdominal_elective_duration_ge_120 = outcome_qc_abdominal,
  baseline_source_qc = baseline_source_qc
)

write_qc_json(qc, "11N_build_aki_outcome.json")

message("==== AKI OUTCOME COMPLETE ====")
message("Saved: ", file.path(DERIVED_DIR, "11N_aki_outcome_creatinine_7d.rds"))