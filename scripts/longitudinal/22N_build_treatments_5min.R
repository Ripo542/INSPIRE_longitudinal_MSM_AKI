# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 22N_build_treatments_5min.R
# Purpose:
#   Build 5-min interval-level treatment exposures from INSPIRE vitals.csv.
#   Treatment records are treated as interval-level events, not continuous signals.
# =============================================================================

SCRIPT <- "22N_build_treatments_5min"

source("R/00_setup/00N_setup_project.R")

message("==== BUILD 5-MIN TREATMENT EXPOSURES ====")

# -----------------------------------------------------------------------------
# LOAD BASE WINDOWS
# -----------------------------------------------------------------------------

base <- readRDS(file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds")) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id)
  ) %>%
  dplyr::filter(has_valid_window) %>%
  dplyr::select(op_id, subject_id, t0, t1)

message("Base valid operations: ", dplyr::n_distinct(base$op_id))

# -----------------------------------------------------------------------------
# LOAD TREATMENT RECORDS
# -----------------------------------------------------------------------------

vitals_path <- file.path(RAW_DIR, "vitals.csv")
stopifnot(file.exists(vitals_path))

fluid_items <- c(
  "ns", "hes", "alb5", "alb20",
  "rbc", "ffp", "cryo"
)

crystalloid_items <- c("ns")
colloid_items <- c("hes", "alb5", "alb20")
blood_items <- c("rbc", "ffp", "cryo")

vasopressor_items <- c("eph", "phe", "nepi")
other_vasoactive_items <- c("epi", "vaso", "dobui", "dopai", "ntgi", "mlni")

treatment_items <- c(
  fluid_items,
  vasopressor_items,
  other_vasoactive_items
)

tx <- data.table::fread(
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

message("Treatment rows before intraop restriction: ", nrow(tx))
message("Operations with treatment rows before restriction: ", dplyr::n_distinct(tx$op_id))

# -----------------------------------------------------------------------------
# RESTRICT TO INTRAOPERATIVE WINDOW AND BIN TO 5-MIN INTERVALS
# -----------------------------------------------------------------------------

tx_intraop <- tx %>%
  dplyr::inner_join(
    base %>% dplyr::select(op_id, t0, t1),
    by = "op_id"
  ) %>%
  dplyr::mutate(
    t_rel = chart_time - as.numeric(t0)
  ) %>%
  dplyr::filter(
    !is.na(t_rel),
    t_rel >= 0,
    chart_time <= as.numeric(t1)
  ) %>%
  dplyr::mutate(
    t = floor(t_rel / 5) * 5,
    value_positive = dplyr::if_else(!is.na(value) & value > 0, 1L, 0L),
    item_group = dplyr::case_when(
      item_name %in% crystalloid_items ~ "crystalloid",
      item_name %in% colloid_items ~ "colloid",
      item_name %in% blood_items ~ "blood",
      item_name %in% vasopressor_items ~ "vasopressor",
      item_name %in% other_vasoactive_items ~ "other_vasoactive",
      TRUE ~ "other"
    )
  )

message("Treatment rows within intraop window: ", nrow(tx_intraop))
message("Operations with intraop treatment rows: ", dplyr::n_distinct(tx_intraop$op_id))

# -----------------------------------------------------------------------------
# INTERVAL-LEVEL AGGREGATION
# -----------------------------------------------------------------------------

tx_interval <- tx_intraop %>%
  dplyr::group_by(op_id, t) %>%
  dplyr::summarise(
    fluid_any = as.integer(any(item_name %in% fluid_items & value_positive == 1L, na.rm = TRUE)),
    crystalloid_any = as.integer(any(item_name %in% crystalloid_items & value_positive == 1L, na.rm = TRUE)),
    colloid_any = as.integer(any(item_name %in% colloid_items & value_positive == 1L, na.rm = TRUE)),
    blood_any = as.integer(any(item_name %in% blood_items & value_positive == 1L, na.rm = TRUE)),
    
    fluid_total_value = sum(dplyr::if_else(item_name %in% fluid_items, value, 0), na.rm = TRUE),
    crystalloid_total_value = sum(dplyr::if_else(item_name %in% crystalloid_items, value, 0), na.rm = TRUE),
    colloid_total_value = sum(dplyr::if_else(item_name %in% colloid_items, value, 0), na.rm = TRUE),
    blood_total_value = sum(dplyr::if_else(item_name %in% blood_items, value, 0), na.rm = TRUE),
    
    vasopressor_any = as.integer(any(item_name %in% vasopressor_items & value_positive == 1L, na.rm = TRUE)),
    eph_any = as.integer(any(item_name == "eph" & value_positive == 1L, na.rm = TRUE)),
    phe_any = as.integer(any(item_name == "phe" & value_positive == 1L, na.rm = TRUE)),
    nepi_any = as.integer(any(item_name == "nepi" & value_positive == 1L, na.rm = TRUE)),
    
    eph_value = sum(dplyr::if_else(item_name == "eph", value, 0), na.rm = TRUE),
    phe_value = sum(dplyr::if_else(item_name == "phe", value, 0), na.rm = TRUE),
    nepi_value = sum(dplyr::if_else(item_name == "nepi", value, 0), na.rm = TRUE),
    
    other_vasoactive_any = as.integer(any(item_name %in% other_vasoactive_items & value_positive == 1L, na.rm = TRUE)),
    
    n_treatment_records = dplyr::n(),
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# COMPLETE AGAINST LONGITUDINAL INTERVAL GRID
# -----------------------------------------------------------------------------

long_base <- readRDS(file.path(DERIVED_DIR, "20N_longitudinal_base_intervals.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

tx_long <- long_base %>%
  dplyr::select(op_id, subject_id, t) %>%
  dplyr::left_join(tx_interval, by = c("op_id", "t")) %>%
  dplyr::mutate(
    dplyr::across(
      c(
        fluid_any, crystalloid_any, colloid_any, blood_any,
        vasopressor_any, eph_any, phe_any, nepi_any,
        other_vasoactive_any
      ),
      ~ dplyr::coalesce(.x, 0L)
    ),
    dplyr::across(
      c(
        fluid_total_value, crystalloid_total_value, colloid_total_value,
        blood_total_value, eph_value, phe_value, nepi_value,
        n_treatment_records
      ),
      ~ dplyr::coalesce(.x, 0)
    )
  ) %>%
  dplyr::arrange(op_id, t) %>%
  dplyr::group_by(op_id) %>%
  dplyr::mutate(
    fluid_any_prev = dplyr::lag(fluid_any, default = 0L),
    vasopressor_any_prev = dplyr::lag(vasopressor_any, default = 0L),
    
    cumulative_fluid_value_prev = dplyr::lag(cumsum(fluid_total_value), default = 0),
    cumulative_vasopressor_intervals_prev = dplyr::lag(cumsum(vasopressor_any), default = 0),
    
    cumulative_crystalloid_value_prev = dplyr::lag(cumsum(crystalloid_total_value), default = 0),
    cumulative_colloid_value_prev = dplyr::lag(cumsum(colloid_total_value), default = 0),
    cumulative_blood_value_prev = dplyr::lag(cumsum(blood_total_value), default = 0)
  ) %>%
  dplyr::ungroup()

# -----------------------------------------------------------------------------
# QC
# -----------------------------------------------------------------------------

qc <- tx_long %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    
    fluid_interval_pct = 100 * mean(fluid_any == 1, na.rm = TRUE),
    crystalloid_interval_pct = 100 * mean(crystalloid_any == 1, na.rm = TRUE),
    colloid_interval_pct = 100 * mean(colloid_any == 1, na.rm = TRUE),
    blood_interval_pct = 100 * mean(blood_any == 1, na.rm = TRUE),
    
    vasopressor_interval_pct = 100 * mean(vasopressor_any == 1, na.rm = TRUE),
    eph_interval_pct = 100 * mean(eph_any == 1, na.rm = TRUE),
    phe_interval_pct = 100 * mean(phe_any == 1, na.rm = TRUE),
    nepi_interval_pct = 100 * mean(nepi_any == 1, na.rm = TRUE),
    
    other_vasoactive_interval_pct = 100 * mean(other_vasoactive_any == 1, na.rm = TRUE)
  )

message("==== Treatment interval QC ====")
print(qc)

op_qc <- tx_long %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    any_fluid = any(fluid_any == 1),
    any_crystalloid = any(crystalloid_any == 1),
    any_colloid = any(colloid_any == 1),
    any_blood = any(blood_any == 1),
    
    any_vasopressor = any(vasopressor_any == 1),
    any_eph = any(eph_any == 1),
    any_phe = any(phe_any == 1),
    any_nepi = any(nepi_any == 1),
    
    total_fluid_value = sum(fluid_total_value, na.rm = TRUE),
    total_crystalloid_value = sum(crystalloid_total_value, na.rm = TRUE),
    total_colloid_value = sum(colloid_total_value, na.rm = TRUE),
    total_blood_value = sum(blood_total_value, na.rm = TRUE),
    
    total_vasopressor_intervals = sum(vasopressor_any, na.rm = TRUE),
    .groups = "drop"
  )

op_summary <- op_qc %>%
  dplyr::summarise(
    n_ops = dplyr::n(),
    
    ops_any_fluid = sum(any_fluid),
    ops_any_crystalloid = sum(any_crystalloid),
    ops_any_colloid = sum(any_colloid),
    ops_any_blood = sum(any_blood),
    
    ops_any_vasopressor = sum(any_vasopressor),
    ops_any_eph = sum(any_eph),
    ops_any_phe = sum(any_phe),
    ops_any_nepi = sum(any_nepi),
    
    total_fluid_value_median = stats::median(total_fluid_value, na.rm = TRUE),
    total_fluid_value_p75 = as.numeric(stats::quantile(total_fluid_value, 0.75, na.rm = TRUE)),
    total_fluid_value_p99 = as.numeric(stats::quantile(total_fluid_value, 0.99, na.rm = TRUE)),
    
    total_vasopressor_intervals_median = stats::median(total_vasopressor_intervals, na.rm = TRUE),
    total_vasopressor_intervals_p75 = as.numeric(stats::quantile(total_vasopressor_intervals, 0.75, na.rm = TRUE)),
    total_vasopressor_intervals_p99 = as.numeric(stats::quantile(total_vasopressor_intervals, 0.99, na.rm = TRUE))
  )

message("==== Treatment operation-level QC ====")
print(op_summary)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(
  tx_long,
  file.path(DERIVED_DIR, "22N_treatments_5min_long.rds")
)

saveRDS(
  op_qc,
  file.path(DERIVED_DIR, "22N_treatments_operation_qc.rds")
)

readr::write_csv(
  qc,
  file.path(QC_DIR, "22N_treatment_interval_qc.csv")
)

readr::write_csv(
  op_summary,
  file.path(QC_DIR, "22N_treatment_operation_qc.csv")
)

qc_json <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  n_rows = qc$n_rows,
  n_ops = qc$n_ops,
  interval_qc = qc,
  operation_qc = op_summary
)

write_qc_json(qc_json, "22N_build_treatments_5min.json")

message("==== 5-MIN TREATMENT EXPOSURES COMPLETE ====")
message("Saved: ", file.path(DERIVED_DIR, "22N_treatments_5min_long.rds"))