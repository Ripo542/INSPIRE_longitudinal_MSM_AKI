# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 20N_build_longitudinal_intervals.R
# Purpose:
#   Build 5-min longitudinal dataset with MAP and HR.
# =============================================================================

SCRIPT <- "20N_build_longitudinal_intervals"

source("R/00_setup/00N_setup_project.R")

message("==== BUILD LONGITUDINAL INTERVALS ====")

# -----------------------------------------------------------------------------
# LOAD BASE COHORT
# -----------------------------------------------------------------------------

base <- readRDS(file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

base_valid <- base %>%
  dplyr::filter(has_valid_window) %>%
  dplyr::mutate(op_id = as.character(op_id)) %>%
  dplyr::select(op_id, subject_id, t0, t1)

message("Base valid operations: ", dplyr::n_distinct(base_valid$op_id))

# -----------------------------------------------------------------------------
# LOAD VITALS
# -----------------------------------------------------------------------------

vitals_path <- file.path(RAW_DIR, "vitals.csv")
stopifnot(file.exists(vitals_path))

vitals <- data.table::fread(
  vitals_path,
  select = c("op_id", "chart_time", "item_name", "value"),
  showProgress = TRUE
) %>%
  janitor::clean_names() %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    chart_time = as.numeric(chart_time),
    value = as.numeric(value)
  )

message("Vitals rows loaded: ", nrow(vitals))
message("op_id class in base_valid: ", paste(class(base_valid$op_id), collapse = ", "))
message("op_id class in vitals: ", paste(class(vitals$op_id), collapse = ", "))

# -----------------------------------------------------------------------------
# CORE SIGNALS
# -----------------------------------------------------------------------------

map_items <- c("art_mbp", "nibp_mbp")
hr_item <- "hr"

vitals_core <- vitals %>%
  dplyr::filter(item_name %in% c(map_items, hr_item)) %>%
  dplyr::inner_join(base_valid, by = "op_id") %>%
  dplyr::mutate(
    t_rel = chart_time - as.numeric(t0)
  ) %>%
  dplyr::filter(
    !is.na(t_rel),
    t_rel >= 0,
    chart_time <= as.numeric(t1)
  )

message("Vitals core rows within intraop window: ", nrow(vitals_core))
message("Operations with core vitals: ", dplyr::n_distinct(vitals_core$op_id))

# -----------------------------------------------------------------------------
# BUILD 5-MIN INTERVAL GRID
# -----------------------------------------------------------------------------

intervals <- base_valid %>%
  dplyr::mutate(
    duration_min = as.numeric(t1 - t0),
    max_bin = floor(duration_min / 5) * 5
  ) %>%
  dplyr::filter(!is.na(max_bin), max_bin >= 0) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(t = list(seq(0, max_bin, by = 5))) %>%
  tidyr::unnest(t) %>%
  dplyr::ungroup() %>%
  dplyr::select(op_id, subject_id, t0, t1, t)

message("Total intervals: ", nrow(intervals))
message("Operations with intervals: ", dplyr::n_distinct(intervals$op_id))

# -----------------------------------------------------------------------------
# BIN VITALS TO 5-MIN INTERVALS
# -----------------------------------------------------------------------------

vitals_core <- vitals_core %>%
  dplyr::mutate(
    t_bin = floor(t_rel / 5) * 5
  )

# MAP: arterial priority over NIBP if both exist in same interval.
map_data <- vitals_core %>%
  dplyr::filter(item_name %in% map_items) %>%
  dplyr::mutate(
    map_priority = dplyr::case_when(
      item_name == "art_mbp" ~ 1L,
      item_name == "nibp_mbp" ~ 2L,
      TRUE ~ 9L
    )
  ) %>%
  dplyr::arrange(op_id, t_bin, map_priority) %>%
  dplyr::group_by(op_id, t_bin) %>%
  dplyr::summarise(
    map = dplyr::first(value),
    map_source = dplyr::first(item_name),
    .groups = "drop"
  )

hr_data <- vitals_core %>%
  dplyr::filter(item_name == hr_item) %>%
  dplyr::group_by(op_id, t_bin) %>%
  dplyr::summarise(
    hr = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# MERGE
# -----------------------------------------------------------------------------

long <- intervals %>%
  dplyr::left_join(map_data, by = c("op_id", "t" = "t_bin")) %>%
  dplyr::left_join(hr_data, by = c("op_id", "t" = "t_bin")) %>%
  dplyr::arrange(op_id, t) %>%
  dplyr::group_by(op_id) %>%
  dplyr::mutate(
    interval_id = dplyr::row_number(),
    time_elapsed_min = t,
    hypotension65 = dplyr::if_else(!is.na(map) & map < 65, 1L, 0L, missing = NA_integer_),
    map_prev = dplyr::lag(map),
    hr_prev = dplyr::lag(hr),
    hypotension65_prev = dplyr::lag(hypotension65),
    map_delta_prev = map - dplyr::lag(map)
  ) %>%
  dplyr::ungroup()

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

qc <- long %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    map_missing_pct = 100 * mean(is.na(map)),
    hr_missing_pct = 100 * mean(is.na(hr)),
    arterial_map_pct = 100 * mean(map_source == "art_mbp", na.rm = TRUE),
    nibp_map_pct = 100 * mean(map_source == "nibp_mbp", na.rm = TRUE),
    hypotension65_pct = 100 * mean(hypotension65 == 1, na.rm = TRUE),
    interval_median_per_op = stats::median(as.numeric(table(op_id)), na.rm = TRUE)
  )

message("==== Longitudinal interval QC ====")
print(qc)

coverage_by_op <- long %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    n_intervals = dplyr::n(),
    map_nonmissing_frac = mean(!is.na(map)),
    hr_nonmissing_frac = mean(!is.na(hr)),
    any_map = any(!is.na(map)),
    any_hr = any(!is.na(hr)),
    .groups = "drop"
  )

coverage_qc <- coverage_by_op %>%
  dplyr::summarise(
    n_ops = dplyr::n(),
    ops_any_map = sum(any_map),
    ops_any_hr = sum(any_hr),
    map_nonmissing_frac_median = stats::median(map_nonmissing_frac, na.rm = TRUE),
    map_nonmissing_frac_p25 = as.numeric(stats::quantile(map_nonmissing_frac, 0.25, na.rm = TRUE)),
    map_nonmissing_frac_p75 = as.numeric(stats::quantile(map_nonmissing_frac, 0.75, na.rm = TRUE)),
    hr_nonmissing_frac_median = stats::median(hr_nonmissing_frac, na.rm = TRUE)
  )

message("==== Coverage by operation QC ====")
print(coverage_qc)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  long,
  file.path(DERIVED_DIR, "20N_longitudinal_base_intervals.rds")
)

saveRDS(
  coverage_by_op,
  file.path(DERIVED_DIR, "20N_longitudinal_coverage_by_op.rds")
)

readr::write_csv(
  qc,
  file.path(QC_DIR, "20N_longitudinal_interval_qc.csv")
)

readr::write_csv(
  coverage_qc,
  file.path(QC_DIR, "20N_longitudinal_coverage_qc.csv")
)

qc_json <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  n_rows = qc$n_rows,
  n_ops = qc$n_ops,
  map_missing_pct = qc$map_missing_pct,
  hr_missing_pct = qc$hr_missing_pct,
  arterial_map_pct = qc$arterial_map_pct,
  nibp_map_pct = qc$nibp_map_pct,
  hypotension65_pct = qc$hypotension65_pct,
  coverage_qc = coverage_qc
)

write_qc_json(qc_json, "20N_build_longitudinal_intervals.json")

message("==== LONGITUDINAL INTERVAL DATASET COMPLETE ====")
message("Saved: ", file.path(DERIVED_DIR, "20N_longitudinal_base_intervals.rds"))