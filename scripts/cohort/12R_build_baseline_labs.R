# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 12R_build_baseline_labs.R
# Purpose:
#   Build baseline haemoglobin and albumin variables for revised MSM pipeline.
#   Baseline definition:
#     1) closest preoperative value within 30 days before t0
#     2) fallback earliest perioperative value from t0 to 6 h after t0
# =============================================================================

SCRIPT <- "12R_build_baseline_labs"

source("R/00_setup/00N_setup_project.R")

message("==== BUILD BASELINE LABS: HAEMOGLOBIN + ALBUMIN ====")

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
    item_name = as.character(item_name),
    value = as.numeric(value)
  ) %>%
  dplyr::filter(!is.na(value))

message("Labs rows loaded: ", nrow(labs))

message("Available lab item names:")
print(sort(unique(labs$item_name)))

# -----------------------------------------------------------------------------
# DEFINE LAB ITEMS
# -----------------------------------------------------------------------------
# Use exact item names after inspection. If item names differ, this script stops.

lab_map <- tibble::tribble(
  ~lab_variable,          ~item_name,
  "baseline_haemoglobin", "hb",
  "baseline_albumin",     "albumin"
)

missing_items <- setdiff(lab_map$item_name, unique(labs$item_name))

if (length(missing_items) > 0) {
  stop(
    "Expected lab item(s) not found in labs.csv: ",
    paste(missing_items, collapse = ", "),
    "\nInspect printed lab item names and update lab_map."
  )
}

labs_sel <- labs %>%
  dplyr::inner_join(lab_map, by = "item_name") %>%
  dplyr::filter(value > 0)

# -----------------------------------------------------------------------------
# QC LAB ITEM SUMMARY
# -----------------------------------------------------------------------------

lab_item_summary <- labs_sel %>%
  dplyr::group_by(lab_variable, item_name) %>%
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

message("==== Baseline lab item summary ====")
print(lab_item_summary)

# -----------------------------------------------------------------------------
# JOIN LABS TO OPERATIONS
# -----------------------------------------------------------------------------

labs_op <- base %>%
  dplyr::select(op_id, subject_id, t0) %>%
  dplyr::inner_join(labs_sel, by = "subject_id") %>%
  dplyr::mutate(
    rel_to_t0_min = chart_time - t0
  )

message("Operation-lab rows after subject join: ", nrow(labs_op))
message("Procedures with any selected lab: ", dplyr::n_distinct(labs_op$op_id))

# -----------------------------------------------------------------------------
# FUNCTION TO BUILD BASELINE LAB
# -----------------------------------------------------------------------------

build_baseline_lab <- function(data, lab_name) {
  
  dat <- data %>%
    dplyr::filter(lab_variable == lab_name)
  
  preop <- dat %>%
    dplyr::filter(rel_to_t0_min < 0, rel_to_t0_min >= -30 * 24 * 60) %>%
    dplyr::mutate(abs_time_to_t0 = abs(rel_to_t0_min)) %>%
    dplyr::arrange(op_id, abs_time_to_t0) %>%
    dplyr::group_by(op_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      op_id,
      value_preop = value,
      time_preop = chart_time,
      source_preop = "preoperative_30d",
      rel_preop = rel_to_t0_min
    )
  
  early <- dat %>%
    dplyr::filter(rel_to_t0_min >= 0, rel_to_t0_min <= 6 * 60) %>%
    dplyr::arrange(op_id, rel_to_t0_min) %>%
    dplyr::group_by(op_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      op_id,
      value_early = value,
      time_early = chart_time,
      source_early = "early_perioperative_0_6h",
      rel_early = rel_to_t0_min
    )
  
  out <- base %>%
    dplyr::select(op_id) %>%
    dplyr::left_join(preop, by = "op_id") %>%
    dplyr::left_join(early, by = "op_id") %>%
    dplyr::mutate(
      value = dplyr::coalesce(value_preop, value_early),
      time = dplyr::coalesce(time_preop, time_early),
      source = dplyr::coalesce(source_preop, source_early),
      rel_to_t0_min = dplyr::coalesce(rel_preop, rel_early)
    ) %>%
    dplyr::transmute(
      op_id,
      !!lab_name := value,
      !!paste0(lab_name, "_time") := time,
      !!paste0(lab_name, "_source") := source,
      !!paste0(lab_name, "_rel_to_t0_min") := rel_to_t0_min
    )
  
  out
}

hb <- build_baseline_lab(labs_op, "baseline_haemoglobin")
alb <- build_baseline_lab(labs_op, "baseline_albumin")

baseline_labs <- base %>%
  dplyr::select(op_id, subject_id) %>%
  dplyr::left_join(hb, by = "op_id") %>%
  dplyr::left_join(alb, by = "op_id") %>%
  dplyr::mutate(
    has_baseline_haemoglobin = !is.na(baseline_haemoglobin),
    has_baseline_albumin = !is.na(baseline_albumin)
  )

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

qc_all <- baseline_labs %>%
  dplyr::summarise(
    n_procedures = dplyr::n(),
    haemoglobin_available = sum(has_baseline_haemoglobin, na.rm = TRUE),
    haemoglobin_available_pct = 100 * mean(has_baseline_haemoglobin, na.rm = TRUE),
    albumin_available = sum(has_baseline_albumin, na.rm = TRUE),
    albumin_available_pct = 100 * mean(has_baseline_albumin, na.rm = TRUE),
    haemoglobin_median = stats::median(baseline_haemoglobin, na.rm = TRUE),
    haemoglobin_p25 = as.numeric(stats::quantile(baseline_haemoglobin, 0.25, na.rm = TRUE)),
    haemoglobin_p75 = as.numeric(stats::quantile(baseline_haemoglobin, 0.75, na.rm = TRUE)),
    albumin_median = stats::median(baseline_albumin, na.rm = TRUE),
    albumin_p25 = as.numeric(stats::quantile(baseline_albumin, 0.25, na.rm = TRUE)),
    albumin_p75 = as.numeric(stats::quantile(baseline_albumin, 0.75, na.rm = TRUE))
  )

hb_source_qc <- baseline_labs %>%
  dplyr::count(baseline_haemoglobin_source, sort = TRUE)

alb_source_qc <- baseline_labs %>%
  dplyr::count(baseline_albumin_source, sort = TRUE)

message("==== Baseline labs QC ====")
print(qc_all)

message("==== Haemoglobin source QC ====")
print(hb_source_qc)

message("==== Albumin source QC ====")
print(alb_source_qc)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  baseline_labs,
  file.path(DERIVED_DIR, "12R_baseline_haemoglobin_albumin.rds")
)

readr::write_csv(
  baseline_labs,
  file.path(TABLES_DIR, "12R_baseline_haemoglobin_albumin.csv")
)

readr::write_csv(
  lab_item_summary,
  file.path(QC_DIR, "12R_baseline_lab_item_summary.csv")
)

readr::write_csv(
  qc_all,
  file.path(QC_DIR, "12R_baseline_labs_qc.csv")
)

readr::write_csv(
  hb_source_qc,
  file.path(QC_DIR, "12R_haemoglobin_source_qc.csv")
)

readr::write_csv(
  alb_source_qc,
  file.path(QC_DIR, "12R_albumin_source_qc.csv")
)

write_qc_json(
  list(
    script = SCRIPT,
    timestamp = as.character(Sys.time()),
    lab_item_summary = lab_item_summary,
    baseline_labs_qc = qc_all,
    haemoglobin_source_qc = hb_source_qc,
    albumin_source_qc = alb_source_qc
  ),
  "12R_build_baseline_labs.json"
)

message("==== BASELINE LABS COMPLETE ====")
message("Saved: ", file.path(DERIVED_DIR, "12R_baseline_haemoglobin_albumin.rds"))