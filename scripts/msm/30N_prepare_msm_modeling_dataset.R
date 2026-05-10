# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 30N_prepare_msm_modeling_dataset.R
# Purpose:
#   Prepare reduced, clean dataset for MSM (treatment models and weights).
#   Select variables, standardise continuous covariates, and check positivity.
# =============================================================================

SCRIPT <- "30N_prepare_msm_modeling_dataset"

source("R/00_setup/00N_setup_project.R")

message("==== PREPARE MSM MODELING DATASET ====")

# -----------------------------------------------------------------------------
# LOAD PRIMARY DATASET
# -----------------------------------------------------------------------------

dat <- readRDS(file.path(FROZEN_DIR, "24N_primary_msm_dataset_duration_ge_120.rds")) %>%
  dplyr::mutate(op_id = as.character(op_id))

message("Rows loaded: ", nrow(dat))
message("Operations: ", dplyr::n_distinct(dat$op_id))

# -----------------------------------------------------------------------------
# SELECT VARIABLES FOR MSM
# -----------------------------------------------------------------------------

msm <- dat %>%
  dplyr::select(
    # identifiers
    op_id, subject_id, interval_id, t,
    
    # outcomes
    aki_creat_7d,
    
    # treatments (current)
    fluid_any,
    vasopressor_any,
    
    # lagged treatments
    fluid_any_lag1,
    vasopressor_any_lag1,
    
    # physiological state (lagged → key confounders)
    map_lag1,
    hr_lag1,
    hypotension65_lag1,
    map_deficit65_lag1,
    map_change_prev,
    
    # cumulative history
    cumulative_hypotension_intervals_prev,
    cumulative_map_deficit65_prev,
    cumulative_fluid_any_intervals_prev,
    cumulative_fluid_value_prev,
    cumulative_vasopressor_intervals_prev,
    
    # time
    time_elapsed_min,
    time_elapsed_hr,
    
    # baseline covariates
    age,
    female,
    asa,
    bmi,
    department_clean,
    baseline_creatinine
  )

# -----------------------------------------------------------------------------
# REMOVE ROWS WITHOUT REQUIRED LAGGED STATE
# -----------------------------------------------------------------------------

msm <- msm %>%
  dplyr::filter(
    !is.na(map_lag1),
    !is.na(hr_lag1),
    !is.na(hypotension65_lag1)
  )

message("Rows after lag filtering: ", nrow(msm))
message("Operations after lag filtering: ", dplyr::n_distinct(msm$op_id))

# -----------------------------------------------------------------------------
# STANDARDISE CONTINUOUS VARIABLES (for model stability)
# -----------------------------------------------------------------------------

scale_vars <- c(
  "map_lag1",
  "hr_lag1",
  "map_deficit65_lag1",
  "map_change_prev",
  "cumulative_hypotension_intervals_prev",
  "cumulative_map_deficit65_prev",
  "cumulative_fluid_value_prev",
  "cumulative_vasopressor_intervals_prev",
  "time_elapsed_hr",
  "age",
  "bmi",
  "baseline_creatinine"
)

scale_safe <- function(x) {
  if (all(is.na(x))) return(x)
  if (sd(x, na.rm = TRUE) == 0) return(x)
  as.numeric(scale(x))
}

msm <- msm %>%
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(scale_vars),
      scale_safe,
      .names = "{.col}_z"
    )
  )

# -----------------------------------------------------------------------------
# FACTOR CODING
# -----------------------------------------------------------------------------

msm <- msm %>%
  dplyr::mutate(
    female = as.integer(female),
    asa = as.factor(asa),
    department_clean = as.factor(department_clean)
  )

# -----------------------------------------------------------------------------
# POSITIVITY CHECK (CRITICAL)
# -----------------------------------------------------------------------------

# Check treatment variability within strata of hypotension

positivity_map <- msm %>%
  dplyr::group_by(hypotension65_lag1) %>%
  dplyr::summarise(
    n = dplyr::n(),
    fluid_rate = mean(fluid_any),
    vasopressor_rate = mean(vasopressor_any),
    .groups = "drop"
  )

message("==== Positivity check: by hypotension ====")
print(positivity_map)

# crude extremes
extreme_fluid <- mean(msm$fluid_any %in% c(0,1))
extreme_vaso  <- mean(msm$vasopressor_any %in% c(0,1))

message("Fluid any rate: ", mean(msm$fluid_any))
message("Vasopressor any rate: ", mean(msm$vasopressor_any))

# -----------------------------------------------------------------------------
# ADD WEIGHT MODEL FORMULAS (for next step)
# -----------------------------------------------------------------------------

# We don’t fit yet, but define structure clearly

formula_fluid <- fluid_any ~
  map_lag1_z +
  hr_lag1_z +
  hypotension65_lag1 +
  map_deficit65_lag1_z +
  map_change_prev_z +
  cumulative_hypotension_intervals_prev_z +
  cumulative_fluid_any_intervals_prev +
  cumulative_vasopressor_intervals_prev_z +
  time_elapsed_hr_z +
  age_z +
  female +
  asa +
  bmi_z +
  baseline_creatinine_z +
  department_clean

formula_vasopressor <- vasopressor_any ~
  map_lag1_z +
  hr_lag1_z +
  hypotension65_lag1 +
  map_deficit65_lag1_z +
  map_change_prev_z +
  cumulative_hypotension_intervals_prev_z +
  cumulative_vasopressor_intervals_prev_z +
  cumulative_fluid_any_intervals_prev +
  time_elapsed_hr_z +
  age_z +
  female +
  asa +
  bmi_z +
  baseline_creatinine_z +
  department_clean

saveRDS(
  list(
    fluid = formula_fluid,
    vasopressor = formula_vasopressor
  ),
  file.path(DERIVED_DIR, "30N_msm_formulas.rds")
)

# -----------------------------------------------------------------------------
# QC SUMMARY
# -----------------------------------------------------------------------------

qc <- msm %>%
  dplyr::summarise(
    n_rows = dplyr::n(),
    n_ops = dplyr::n_distinct(op_id),
    fluid_rate = mean(fluid_any),
    vasopressor_rate = mean(vasopressor_any),
    hypotension_rate = mean(hypotension65_lag1)
  )

message("==== MSM modeling dataset QC ====")
print(qc)

saveRDS(
  msm,
  file.path(DERIVED_DIR, "30N_msm_modeling_dataset.rds")
)

write_qc_json(qc, "30N_prepare_msm_modeling_dataset.json")

message("==== MSM MODELING DATASET READY ====")
message("Saved: ", file.path(DERIVED_DIR, "30N_msm_modeling_dataset.rds"))
