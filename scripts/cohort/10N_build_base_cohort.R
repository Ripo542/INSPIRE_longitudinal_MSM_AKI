# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 10N_build_base_cohort.R
# Purpose:
#   Build base surgical cohort from operations.csv for longitudinal MSM pipeline.
#   Uses INSPIRE v1.3 column names.
# =============================================================================

SCRIPT <- "10N_build_base_cohort"

source("R/00_setup/00N_setup_project.R")

message("==== BUILD BASE COHORT ====")

ops_path <- file.path(RAW_DIR, "operations.csv")
stopifnot(file.exists(ops_path))

ops <- data.table::fread(ops_path, showProgress = TRUE) %>%
  janitor::clean_names()

message("operations rows: ", nrow(ops))
message("operations columns: ", paste(names(ops), collapse = ", "))

# ---- Save column inventory ----
ops_columns <- tibble::tibble(
  column = names(ops),
  class = vapply(ops, function(x) paste(class(x), collapse = ","), character(1)),
  n_missing = vapply(ops, function(x) sum(is.na(x)), integer(1)),
  missing_pct = vapply(ops, function(x) 100 * mean(is.na(x)), numeric(1))
)

readr::write_csv(
  ops_columns,
  file.path(QC_DIR, "10N_operations_column_inventory.csv")
)

# ---- Required IDs ----
required_ids <- c("op_id", "subject_id")
missing_ids <- setdiff(required_ids, names(ops))

if (length(missing_ids) > 0) {
  stop("Missing required ID columns in operations.csv: ", paste(missing_ids, collapse = ", "))
}

# ---- Helper: first existing column ----
first_existing <- function(candidates, data_names = names(ops)) {
  out <- intersect(candidates, data_names)
  if (length(out) == 0) NA_character_ else out[1]
}

# ---- Detect variables ----
age_col <- first_existing(c("age", "age_yr", "age_years"))
sex_col <- first_existing(c("sex", "gender", "female"))
asa_col <- first_existing(c("asa", "asa_ps", "asa_class"))
weight_col <- first_existing(c("weight", "body_weight", "weight_kg"))
height_col <- first_existing(c("height", "body_height", "height_cm"))
bmi_col <- first_existing(c("bmi", "body_mass_index"))
dept_col <- first_existing(c("department", "dept", "surgery_department"))
emop_col <- first_existing(c("emop", "emergency", "emergency_op", "urgent", "urgency"))

# INSPIRE v1.3 time columns
anstart_col <- first_existing(c("anstart_time", "anstart", "anes_start", "anesthesia_start"))
anend_col   <- first_existing(c("anend_time", "anend", "anes_end", "anesthesia_end"))
opstart_col <- first_existing(c("opstart_time", "opstart", "operation_start", "surgery_start"))
opend_col   <- first_existing(c("opend_time", "opend", "operation_end", "surgery_end"))
orin_col    <- first_existing(c("orin_time", "orin", "or_in", "room_in"))
orout_col   <- first_existing(c("orout_time", "orout", "or_out", "room_out"))

detected_cols <- tibble::tibble(
  concept = c(
    "age", "sex", "asa", "weight", "height", "bmi",
    "department", "emergency",
    "anstart", "anend", "opstart", "opend", "orin", "orout"
  ),
  column = c(
    age_col, sex_col, asa_col, weight_col, height_col, bmi_col,
    dept_col, emop_col,
    anstart_col, anend_col, opstart_col, opend_col, orin_col, orout_col
  )
)

message("==== Detected columns ====")
print(detected_cols)

readr::write_csv(
  detected_cols,
  file.path(QC_DIR, "10N_operations_detected_columns.csv")
)

# ---- Construct base dataset ----
base <- ops %>%
  dplyr::transmute(
    op_id = as.character(.data$op_id),
    subject_id = as.character(.data$subject_id),
    
    age = if (!is.na(age_col)) .data[[age_col]] else NA_real_,
    sex_raw = if (!is.na(sex_col)) as.character(.data[[sex_col]]) else NA_character_,
    asa = if (!is.na(asa_col)) .data[[asa_col]] else NA,
    weight = if (!is.na(weight_col)) .data[[weight_col]] else NA_real_,
    height = if (!is.na(height_col)) .data[[height_col]] else NA_real_,
    bmi_raw = if (!is.na(bmi_col)) .data[[bmi_col]] else NA_real_,
    department = if (!is.na(dept_col)) as.character(.data[[dept_col]]) else NA_character_,
    emop = if (!is.na(emop_col)) .data[[emop_col]] else NA,
    
    anstart_time = if (!is.na(anstart_col)) .data[[anstart_col]] else NA_real_,
    anend_time   = if (!is.na(anend_col)) .data[[anend_col]] else NA_real_,
    opstart_time = if (!is.na(opstart_col)) .data[[opstart_col]] else NA_real_,
    opend_time   = if (!is.na(opend_col)) .data[[opend_col]] else NA_real_,
    orin_time    = if (!is.na(orin_col)) .data[[orin_col]] else NA_real_,
    orout_time   = if (!is.na(orout_col)) .data[[orout_col]] else NA_real_
  )

# ---- BMI derivation ----
# INSPIRE height is expected in cm and weight in kg. If height appears in m,
# this block handles both.
base <- base %>%
  dplyr::mutate(
    height_m = dplyr::case_when(
      is.na(height) ~ NA_real_,
      height > 3 ~ height / 100,
      height > 0 & height <= 3 ~ height,
      TRUE ~ NA_real_
    ),
    bmi_derived = dplyr::if_else(
      !is.na(weight) & !is.na(height_m) & height_m > 0,
      weight / (height_m^2),
      NA_real_
    ),
    bmi = dplyr::coalesce(as.numeric(bmi_raw), as.numeric(bmi_derived))
  )

# ---- Harmonise sex/female ----
base <- base %>%
  dplyr::mutate(
    female = dplyr::case_when(
      sex_raw %in% c("F", "Female", "female", "f", "2") ~ 1L,
      sex_raw %in% c("M", "Male", "male", "m", "1") ~ 0L,
      sex_raw %in% c("0") ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# ---- Define intraoperative window ----
# For the longitudinal MSM, anesthesia window is primary because treatment and vitals
# are recorded during anesthetic care. Operation and OR windows are retained for QC.
base <- base %>%
  dplyr::mutate(
    t0 = dplyr::case_when(
      !is.na(anstart_time) & !is.na(anend_time) ~ anstart_time,
      !is.na(opstart_time) & !is.na(opend_time) ~ opstart_time,
      !is.na(orin_time) & !is.na(orout_time) ~ orin_time,
      TRUE ~ NA_real_
    ),
    t1 = dplyr::case_when(
      !is.na(anstart_time) & !is.na(anend_time) ~ anend_time,
      !is.na(opstart_time) & !is.na(opend_time) ~ opend_time,
      !is.na(orin_time) & !is.na(orout_time) ~ orout_time,
      TRUE ~ NA_real_
    ),
    time_window_source = dplyr::case_when(
      !is.na(anstart_time) & !is.na(anend_time) ~ "anesthesia",
      !is.na(opstart_time) & !is.na(opend_time) ~ "operation",
      !is.na(orin_time) & !is.na(orout_time) ~ "or_room",
      TRUE ~ "missing"
    ),
    intraop_duration_min = as.numeric(t1 - t0),
    operation_duration_min = as.numeric(opend_time - opstart_time),
    anesthesia_duration_min = as.numeric(anend_time - anstart_time),
    or_duration_min = as.numeric(orout_time - orin_time)
  )

# ---- Eligibility flags ----
base <- base %>%
  dplyr::mutate(
    adult = !is.na(age) & age >= 18,
    
    elective = dplyr::case_when(
      is.na(emop) ~ NA,
      emop %in% c(0, "0", "N", "No", "no", FALSE) ~ TRUE,
      emop %in% c(1, "1", "Y", "Yes", "yes", TRUE) ~ FALSE,
      TRUE ~ NA
    ),
    
    has_valid_window = !is.na(t0) & !is.na(t1) & intraop_duration_min > 0,
    duration_ge_120 = has_valid_window & intraop_duration_min >= 120,
    asa_vi = as.character(asa) %in% c("6", "VI", "ASA VI"),
    
    department_clean = toupper(trimws(department)),
    
    abdominal_department = dplyr::case_when(
      department_clean %in% c("GS", "OG", "UR") ~ TRUE,
      is.na(department_clean) ~ NA,
      TRUE ~ FALSE
    )
  )

# ---- Main base cohort for downstream longitudinal construction ----
# This is deliberately broad. Final analysis cohort will additionally require:
# valid MAP coverage, AKI ascertainment, baseline creatinine, and covariates.
base_cohort <- base %>%
  dplyr::filter(adult) %>%
  dplyr::filter(has_valid_window) %>%
  dplyr::filter(!asa_vi | is.na(asa_vi))

# ---- More restrictive abdominal elective cohort candidate ----
abdominal_elective_candidate <- base_cohort %>%
  dplyr::filter(abdominal_department) %>%
  dplyr::filter(elective | is.na(elective)) %>%
  dplyr::filter(duration_ge_120)

# ---- QC summaries ----
department_summary <- base %>%
  dplyr::count(department_clean, sort = TRUE)

eligibility_qc <- tibble::tibble(
  step = c(
    "operations_total",
    "adult",
    "valid_time_window",
    "adult_valid_window",
    "adult_valid_window_non_asa_vi",
    "adult_valid_window_duration_ge_120",
    "adult_valid_window_elective",
    "adult_valid_window_abdominal_department",
    "candidate_abdominal_elective_duration_ge_120"
  ),
  n = c(
    nrow(base),
    sum(base$adult, na.rm = TRUE),
    sum(base$has_valid_window, na.rm = TRUE),
    sum(base$adult & base$has_valid_window, na.rm = TRUE),
    nrow(base_cohort),
    sum(base$adult & base$has_valid_window & base$duration_ge_120, na.rm = TRUE),
    sum(base$adult & base$has_valid_window & base$elective, na.rm = TRUE),
    sum(base$adult & base$has_valid_window & base$abdominal_department, na.rm = TRUE),
    nrow(abdominal_elective_candidate)
  )
)

window_qc <- base %>%
  dplyr::filter(has_valid_window) %>%
  dplyr::group_by(time_window_source) %>%
  dplyr::summarise(
    n = dplyr::n(),
    duration_median = stats::median(intraop_duration_min, na.rm = TRUE),
    duration_p25 = as.numeric(stats::quantile(intraop_duration_min, 0.25, na.rm = TRUE)),
    duration_p75 = as.numeric(stats::quantile(intraop_duration_min, 0.75, na.rm = TRUE)),
    duration_p99 = as.numeric(stats::quantile(intraop_duration_min, 0.99, na.rm = TRUE)),
    .groups = "drop"
  )

bmi_qc <- base %>%
  dplyr::summarise(
    n = dplyr::n(),
    weight_missing_pct = 100 * mean(is.na(weight)),
    height_missing_pct = 100 * mean(is.na(height)),
    bmi_missing_pct = 100 * mean(is.na(bmi)),
    bmi_median = stats::median(bmi, na.rm = TRUE),
    bmi_p01 = as.numeric(stats::quantile(bmi, 0.01, na.rm = TRUE)),
    bmi_p99 = as.numeric(stats::quantile(bmi, 0.99, na.rm = TRUE))
  )

message("==== Eligibility QC ====")
print(eligibility_qc)

message("==== Time-window QC ====")
print(window_qc)

message("==== BMI QC ====")
print(bmi_qc)

message("==== Department summary ====")
print(head(department_summary, 30))

# ---- Save outputs ----
saveRDS(base, file.path(DERIVED_DIR, "10N_operations_base_all.rds"))
saveRDS(base_cohort, file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds"))
saveRDS(
  abdominal_elective_candidate,
  file.path(FROZEN_DIR, "10N_candidate_abdominal_elective_duration_ge_120.rds")
)

readr::write_csv(eligibility_qc, file.path(QC_DIR, "10N_base_cohort_eligibility_qc.csv"))
readr::write_csv(window_qc, file.path(QC_DIR, "10N_time_window_qc.csv"))
readr::write_csv(bmi_qc, file.path(QC_DIR, "10N_bmi_qc.csv"))
readr::write_csv(department_summary, file.path(QC_DIR, "10N_department_summary.csv"))

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  raw_dir = RAW_DIR,
  n_operations = nrow(base),
  n_base_cohort = nrow(base_cohort),
  n_candidate_abdominal_elective_duration_ge_120 = nrow(abdominal_elective_candidate),
  detected_columns = detected_cols,
  eligibility_qc = eligibility_qc,
  window_qc = window_qc,
  bmi_qc = bmi_qc
)

write_qc_json(qc, "10N_build_base_cohort.json")

message("==== BASE COHORT COMPLETE ====")
message("Saved base cohort: ", file.path(FROZEN_DIR, "10N_base_cohort_adult_valid_window.rds"))
message("Saved candidate cohort: ", file.path(FROZEN_DIR, "10N_candidate_abdominal_elective_duration_ge_120.rds"))