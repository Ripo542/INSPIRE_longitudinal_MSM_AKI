# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 21N_audit_treatment_items.R
# Purpose:
#   Audit treatment-related items in vitals.csv before constructing longitudinal
#   fluid and vasopressor exposures.
# =============================================================================

SCRIPT <- "21N_audit_treatment_items"

source("R/00_setup/00N_setup_project.R")

message("==== AUDIT TREATMENT ITEMS ====")

vitals_path <- file.path(RAW_DIR, "vitals.csv")
stopifnot(file.exists(vitals_path))

treatment_items <- c(
  # fluids / blood products
  "ns", "hes", "alb5", "alb20", "rbc", "ffp", "cryo",
  "d5w", "d10w", "d50w", "hns", "ds",
  
  # vasopressors / vasoactive drugs
  "eph", "phe", "nepi", "epi", "vaso",
  "dobui", "dopai", "ntgi", "mlni"
)

vitals_tx <- data.table::fread(
  vitals_path,
  select = c("op_id", "subject_id", "chart_time", "item_name", "value"),
  showProgress = TRUE
) %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id),
    chart_time = as.numeric(chart_time),
    value = as.numeric(value)
  ) %>%
  dplyr::filter(item_name %in% treatment_items)

message("Treatment rows loaded: ", nrow(vitals_tx))
message("Operations with any treatment item: ", dplyr::n_distinct(vitals_tx$op_id))

# ---- Overall treatment item summary ----
item_summary <- vitals_tx %>%
  dplyr::group_by(item_name) %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    n_subjects = dplyr::n_distinct(subject_id),
    value_missing_pct = 100 * mean(is.na(value)),
    value_min = suppressWarnings(min(value, na.rm = TRUE)),
    value_p01 = suppressWarnings(as.numeric(stats::quantile(value, 0.01, na.rm = TRUE))),
    value_p25 = suppressWarnings(as.numeric(stats::quantile(value, 0.25, na.rm = TRUE))),
    value_median = suppressWarnings(stats::median(value, na.rm = TRUE)),
    value_p75 = suppressWarnings(as.numeric(stats::quantile(value, 0.75, na.rm = TRUE))),
    value_p99 = suppressWarnings(as.numeric(stats::quantile(value, 0.99, na.rm = TRUE))),
    value_max = suppressWarnings(max(value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_rows))

message("==== Treatment item summary ====")
print(item_summary)

readr::write_csv(
  item_summary,
  file.path(QC_DIR, "21N_treatment_item_summary.csv")
)

# ---- Temporal gap summary for treatment records ----
tx_dt <- vitals_tx %>%
  dplyr::arrange(op_id, item_name, chart_time) %>%
  dplyr::group_by(op_id, item_name) %>%
  dplyr::mutate(dt = chart_time - dplyr::lag(chart_time)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(dt), dt > 0)

tx_dt_summary <- tx_dt %>%
  dplyr::group_by(item_name) %>%
  dplyr::summarise(
    n_gaps = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    dt_min = min(dt, na.rm = TRUE),
    dt_median = stats::median(dt, na.rm = TRUE),
    dt_p75 = as.numeric(stats::quantile(dt, 0.75, na.rm = TRUE)),
    dt_p95 = as.numeric(stats::quantile(dt, 0.95, na.rm = TRUE)),
    dt_p99 = as.numeric(stats::quantile(dt, 0.99, na.rm = TRUE)),
    dt_max = max(dt, na.rm = TRUE),
    pct_gap_le_5 = 100 * mean(dt <= 5, na.rm = TRUE),
    pct_gap_gt_5 = 100 * mean(dt > 5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(item_name)

message("==== Treatment temporal gap summary ====")
print(tx_dt_summary)

readr::write_csv(
  tx_dt_summary,
  file.path(QC_DIR, "21N_treatment_temporal_gap_summary.csv")
)

# ---- Values per operation: crude exposure prevalence ----
tx_by_op <- vitals_tx %>%
  dplyr::group_by(op_id, item_name) %>%
  dplyr::summarise(
    n_records = dplyr::n(),
    total_value = sum(value, na.rm = TRUE),
    any_positive_value = any(value > 0, na.rm = TRUE),
    .groups = "drop"
  )

tx_prevalence <- tx_by_op %>%
  dplyr::group_by(item_name) %>%
  dplyr::summarise(
    n_ops_with_record = dplyr::n_distinct(op_id),
    n_ops_positive = sum(any_positive_value, na.rm = TRUE),
    total_value_median_among_recorded = stats::median(total_value, na.rm = TRUE),
    total_value_p75_among_recorded = as.numeric(stats::quantile(total_value, 0.75, na.rm = TRUE)),
    total_value_p99_among_recorded = as.numeric(stats::quantile(total_value, 0.99, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_ops_with_record))

message("==== Crude treatment prevalence by operation ====")
print(tx_prevalence)

readr::write_csv(
  tx_prevalence,
  file.path(QC_DIR, "21N_treatment_prevalence_by_operation.csv")
)

# ---- Save audit object ----
audit_object <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  raw_dir = RAW_DIR,
  treatment_items = treatment_items,
  item_summary = item_summary,
  tx_dt_summary = tx_dt_summary,
  tx_prevalence = tx_prevalence
)

saveRDS(
  audit_object,
  file.path(DERIVED_DIR, "21N_treatment_item_audit.rds")
)

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  raw_dir = RAW_DIR,
  n_treatment_rows = nrow(vitals_tx),
  n_ops_any_treatment_item = dplyr::n_distinct(vitals_tx$op_id),
  treatment_items_found = sort(unique(vitals_tx$item_name))
)

write_qc_json(qc, "21N_audit_treatment_items.json")

message("==== TREATMENT AUDIT COMPLETE ====")