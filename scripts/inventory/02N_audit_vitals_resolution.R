# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 02N_audit_vitals_resolution.R
# Purpose:
#   Audit INSPIRE vitals.csv temporal structure, item names, and effective
#   resolution of arterial pressure and heart-rate signals.
# =============================================================================

SCRIPT <- "02N_audit_vitals_resolution"

source("R/00_setup/00N_setup_project.R")

message("==== VITALS RESOLUTION AUDIT ====")

vitals_path <- file.path(RAW_DIR, "vitals.csv")
stopifnot(file.exists(vitals_path))

# ---- Read required columns ----
vitals <- data.table::fread(
  vitals_path,
  select = c("op_id", "subject_id", "chart_time", "item_name", "value"),
  showProgress = TRUE
)

vitals <- janitor::clean_names(vitals)

message("Rows: ", nrow(vitals))
message("Unique operations: ", data.table::uniqueN(vitals$op_id))
message("Unique subjects: ", data.table::uniqueN(vitals$subject_id))

# ---- Item inventory ----
item_inventory <- vitals %>%
  dplyr::count(item_name, sort = TRUE) %>%
  dplyr::mutate(percent = 100 * n / sum(n))

readr::write_csv(
  item_inventory,
  file.path(QC_DIR, "02N_vitals_item_inventory.csv")
)

message("==== Item inventory saved ====")

# ---- Candidate hemodynamic items ----
candidate_patterns <- c(
  "mbp", "map", "sbp", "dbp", "hr", "heart", "art", "nibp", "bp"
)

candidate_regex <- paste(candidate_patterns, collapse = "|")

candidate_items <- item_inventory %>%
  dplyr::filter(grepl(candidate_regex, item_name, ignore.case = TRUE))

message("==== Candidate hemodynamic item names ====")
print(candidate_items)

readr::write_csv(
  candidate_items,
  file.path(QC_DIR, "02N_candidate_hemodynamic_items.csv")
)

# ---- Core signals ----
core_items <- c(
  "art_mbp", "art_sbp", "art_dbp",
  "nibp_mbp", "nibp_sbp", "nibp_dbp",
  "hr"
)

available_core_items <- intersect(core_items, unique(vitals$item_name))

message("Available core items: ", paste(available_core_items, collapse = ", "))

if (length(available_core_items) == 0) {
  stop("No expected core hemodynamic items found.")
}

# ---- Signal ranges ----
signal_counts <- vitals %>%
  dplyr::filter(item_name %in% available_core_items) %>%
  dplyr::group_by(item_name) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    n_subjects = dplyr::n_distinct(subject_id),
    value_missing_pct = 100 * mean(is.na(value)),
    value_min = suppressWarnings(min(value, na.rm = TRUE)),
    value_p01 = suppressWarnings(as.numeric(stats::quantile(value, 0.01, na.rm = TRUE))),
    value_median = suppressWarnings(stats::median(value, na.rm = TRUE)),
    value_p99 = suppressWarnings(as.numeric(stats::quantile(value, 0.99, na.rm = TRUE))),
    value_max = suppressWarnings(max(value, na.rm = TRUE)),
    .groups = "drop"
  )

message("==== Core signal counts and ranges ====")
print(signal_counts)

readr::write_csv(
  signal_counts,
  file.path(QC_DIR, "02N_core_signal_counts_and_ranges.csv")
)

# ---- Temporal resolution audit ----
signal_dt <- vitals %>%
  dplyr::filter(item_name %in% available_core_items) %>%
  dplyr::filter(!is.na(chart_time)) %>%
  dplyr::arrange(op_id, item_name, chart_time) %>%
  dplyr::group_by(op_id, item_name) %>%
  dplyr::mutate(
    dt = as.numeric(chart_time - dplyr::lag(chart_time))
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(dt), dt > 0)

dt_summary <- signal_dt %>%
  dplyr::group_by(item_name) %>%
  dplyr::summarise(
    n_gaps = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    dt_min = min(dt, na.rm = TRUE),
    dt_p01 = as.numeric(stats::quantile(dt, 0.01, na.rm = TRUE)),
    dt_p05 = as.numeric(stats::quantile(dt, 0.05, na.rm = TRUE)),
    dt_median = stats::median(dt, na.rm = TRUE),
    dt_p75 = as.numeric(stats::quantile(dt, 0.75, na.rm = TRUE)),
    dt_p95 = as.numeric(stats::quantile(dt, 0.95, na.rm = TRUE)),
    dt_p99 = as.numeric(stats::quantile(dt, 0.99, na.rm = TRUE)),
    dt_max = max(dt, na.rm = TRUE),
    pct_gap_1 = 100 * mean(dt == 1, na.rm = TRUE),
    pct_gap_5 = 100 * mean(dt == 5, na.rm = TRUE),
    pct_gap_le_1 = 100 * mean(dt <= 1, na.rm = TRUE),
    pct_gap_le_5 = 100 * mean(dt <= 5, na.rm = TRUE),
    pct_gap_gt_5 = 100 * mean(dt > 5, na.rm = TRUE),
    .groups = "drop"
  )

message("==== Temporal gap summary ====")
print(dt_summary)

readr::write_csv(
  dt_summary,
  file.path(QC_DIR, "02N_vitals_temporal_gap_summary.csv")
)

# ---- Per-operation coverage ----
coverage_by_op <- vitals %>%
  dplyr::filter(item_name %in% available_core_items) %>%
  dplyr::filter(!is.na(chart_time)) %>%
  dplyr::group_by(op_id, item_name) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    first_time = min(chart_time, na.rm = TRUE),
    last_time = max(chart_time, na.rm = TRUE),
    observed_span = as.numeric(last_time - first_time),
    median_gap = {
      x <- sort(unique(chart_time))
      if (length(x) >= 2) stats::median(diff(x), na.rm = TRUE) else NA_real_
    },
    .groups = "drop"
  )

coverage_summary <- coverage_by_op %>%
  dplyr::group_by(item_name) %>%
  dplyr::summarise(
    n_ops = dplyr::n_distinct(op_id),
    n_obs_median = stats::median(n_obs, na.rm = TRUE),
    n_obs_p25 = as.numeric(stats::quantile(n_obs, 0.25, na.rm = TRUE)),
    n_obs_p75 = as.numeric(stats::quantile(n_obs, 0.75, na.rm = TRUE)),
    span_median = stats::median(observed_span, na.rm = TRUE),
    median_gap_median = stats::median(median_gap, na.rm = TRUE),
    median_gap_p25 = as.numeric(stats::quantile(median_gap, 0.25, na.rm = TRUE)),
    median_gap_p75 = as.numeric(stats::quantile(median_gap, 0.75, na.rm = TRUE)),
    .groups = "drop"
  )

message("==== Per-operation coverage summary ====")
print(coverage_summary)

readr::write_csv(
  coverage_summary,
  file.path(QC_DIR, "02N_vitals_per_operation_coverage_summary.csv")
)

# ---- MAP resolution decision ----
map_items <- intersect(c("art_mbp", "nibp_mbp"), available_core_items)

map_dt_summary <- dt_summary %>%
  dplyr::filter(item_name %in% map_items)

map_resolution_decision <- dplyr::case_when(
  nrow(map_dt_summary) == 0 ~
    "No MAP signal found.",
  all(map_dt_summary$dt_median <= 1, na.rm = TRUE) ~
    "MAP appears compatible with minute-level analysis, but 5-min intervals remain primary for reviewer-facing robustness.",
  all(map_dt_summary$dt_median <= 5, na.rm = TRUE) ~
    "MAP is compatible with 5-min interval analysis. One-minute analysis should not be primary.",
  TRUE ~
    "MAP appears coarser than 5-min for at least one source; longitudinal analysis may require coarser intervals or stricter coverage rules."
)

message("==== Resolution decision ====")
message(map_resolution_decision)

# ---- Save audit object ----
audit_object <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  raw_dir = RAW_DIR,
  n_rows = nrow(vitals),
  n_ops = data.table::uniqueN(vitals$op_id),
  n_subjects = data.table::uniqueN(vitals$subject_id),
  item_inventory = item_inventory,
  candidate_items = candidate_items,
  available_core_items = available_core_items,
  signal_counts = signal_counts,
  dt_summary = dt_summary,
  coverage_summary = coverage_summary,
  map_resolution_decision = map_resolution_decision
)

saveRDS(
  audit_object,
  file.path(DERIVED_DIR, "02N_vitals_resolution_audit.rds")
)

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  raw_dir = RAW_DIR,
  n_rows = nrow(vitals),
  n_ops = data.table::uniqueN(vitals$op_id),
  n_subjects = data.table::uniqueN(vitals$subject_id),
  available_core_items = available_core_items,
  map_resolution_decision = map_resolution_decision
)

write_qc_json(qc, "02N_audit_vitals_resolution.json")

message("==== VITALS RESOLUTION AUDIT COMPLETE ====")