# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 40R_build_table2_absolute_risk_differences_revised.R
# Purpose:
#   Build revised Table 2: absolute risk differences for:
#     1) primary elective abdominal major cohort
#     2) sensitivity full major noncardiac cohort
#   using cluster bootstrap with model refitting.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== BUILD REVISED TABLE 2: ABSOLUTE RISK DIFFERENCES ====")

# -----------------------------------------------------------------------------
# Load datasets and models
# -----------------------------------------------------------------------------

dat_abd <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_elective_abdominal_major.rds"))
dat_full <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_full_major_noncardiac.rds"))

fit_abd <- readRDS(file.path(DERIVED_DIR, "34R_model_primary_abdominal.rds"))
fit_full <- readRDS(file.path(DERIVED_DIR, "34R_model_sensitivity_full.rds"))

# -----------------------------------------------------------------------------
# Model formula
# -----------------------------------------------------------------------------

model_formula <- aki_creat_7d ~
  hypotension_time_10min +
  map_deficit65_per_10 +
  fluid_intervals_per_5 +
  vasopressor_intervals_per_5 +
  age_z +
  female +
  asa_final +
  bmi_z +
  baseline_creatinine_z +
  baseline_haemoglobin_z +
  baseline_albumin_z +
  department_final +
  n_intervals_z +
  diabetes +
  ckd +
  ischaemic_heart_disease +
  chronic_pulmonary_disease

# -----------------------------------------------------------------------------
# Helper
# -----------------------------------------------------------------------------

prepare_dat <- function(dat) {
  dat %>%
    dplyr::mutate(
      asa_final = as.factor(asa_final),
      department_final = as.factor(department_final),
      diabetes = as.integer(diabetes),
      ckd = as.integer(ckd),
      ischaemic_heart_disease = as.integer(ischaemic_heart_disease),
      chronic_pulmonary_disease = as.integer(chronic_pulmonary_disease)
    )
}

dat_abd <- prepare_dat(dat_abd)
dat_full <- prepare_dat(dat_full)

compute_ard_boot <- function(dat, fit, cohort_label, var, exposure_label, n_boot = 1000, seed = 20260505) {
  
  message("Computing ARD: ", cohort_label, " | ", exposure_label)
  
  p25 <- as.numeric(stats::quantile(dat[[var]], 0.25, na.rm = TRUE))
  p75 <- as.numeric(stats::quantile(dat[[var]], 0.75, na.rm = TRUE))
  
  low_dat <- dat
  high_dat <- dat
  
  low_dat[[var]] <- p25
  high_dat[[var]] <- p75
  
  p_low <- predict(fit, newdata = low_dat, type = "response")
  p_high <- predict(fit, newdata = high_dat, type = "response")
  
  ard <- mean(p_high - p_low, na.rm = TRUE)
  
  set.seed(seed)
  
  boot <- rep(NA_real_, n_boot)
  
  for (b in seq_len(n_boot)) {
    
    idx <- sample(seq_len(nrow(dat)), replace = TRUE)
    boot_dat <- dat[idx, , drop = FALSE]
    
    boot_dat$asa_final <- factor(boot_dat$asa_final, levels = levels(dat$asa_final))
    boot_dat$department_final <- factor(boot_dat$department_final, levels = levels(dat$department_final))
    
    fit_boot <- tryCatch(
      glm(
        model_formula,
        data = boot_dat,
        weights = sw_1_99,
        family = quasibinomial()
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit_boot)) next
    
    low_boot <- boot_dat
    high_boot <- boot_dat
    
    low_boot[[var]] <- p25
    high_boot[[var]] <- p75
    
    p_low_boot <- tryCatch(
      predict(fit_boot, newdata = low_boot, type = "response"),
      error = function(e) rep(NA_real_, nrow(low_boot))
    )
    
    p_high_boot <- tryCatch(
      predict(fit_boot, newdata = high_boot, type = "response"),
      error = function(e) rep(NA_real_, nrow(high_boot))
    )
    
    boot[b] <- mean(p_high_boot - p_low_boot, na.rm = TRUE)
  }
  
  boot <- boot[is.finite(boot)]
  
  tibble::tibble(
    Cohort = cohort_label,
    Exposure = exposure_label,
    Variable = var,
    Contrast = "25th to 75th percentile",
    Low_value = p25,
    High_value = p75,
    Predicted_risk_low = mean(p_low, na.rm = TRUE),
    Predicted_risk_high = mean(p_high, na.rm = TRUE),
    Absolute_risk_difference = ard,
    CI_low = as.numeric(stats::quantile(boot, 0.025, na.rm = TRUE)),
    CI_high = as.numeric(stats::quantile(boot, 0.975, na.rm = TRUE)),
    Bootstrap_iterations = length(boot)
  )
}

# -----------------------------------------------------------------------------
# Build table
# -----------------------------------------------------------------------------

exposures <- list(
  list(var = "hypotension_time_10min", exposure = "Hypotension duration"),
  list(var = "map_deficit65_per_10", exposure = "MAP deficit <65 mmHg"),
  list(var = "fluid_intervals_per_5", exposure = "Fluid exposure"),
  list(var = "vasopressor_intervals_per_5", exposure = "Vasopressor exposure")
)

table2 <- dplyr::bind_rows(
  lapply(seq_along(exposures), function(i) {
    compute_ard_boot(
      dat = dat_abd,
      fit = fit_abd,
      cohort_label = "Primary: elective major abdominal surgery",
      var = exposures[[i]]$var,
      exposure_label = exposures[[i]]$exposure,
      n_boot = 1000,
      seed = 20260500 + i
    )
  }),
  lapply(seq_along(exposures), function(i) {
    compute_ard_boot(
      dat = dat_full,
      fit = fit_full,
      cohort_label = "Sensitivity: full major noncardiac surgery",
      var = exposures[[i]]$var,
      exposure_label = exposures[[i]]$exposure,
      n_boot = 1000,
      seed = 20260600 + i
    )
  })
)

table2_formatted <- table2 %>%
  dplyr::mutate(
    Low_value = round(Low_value, 2),
    High_value = round(High_value, 2),
    Predicted_risk_low = sprintf("%.1f%%", 100 * Predicted_risk_low),
    Predicted_risk_high = sprintf("%.1f%%", 100 * Predicted_risk_high),
    Absolute_risk_difference = sprintf(
      "%+.2f%% (%.2f to %.2f)",
      100 * Absolute_risk_difference,
      100 * CI_low,
      100 * CI_high
    )
  ) %>%
  dplyr::select(
    Cohort,
    Exposure,
    Contrast,
    Low_value,
    High_value,
    Predicted_risk_low,
    Predicted_risk_high,
    Absolute_risk_difference,
    Bootstrap_iterations
  )

message("==== Revised Table 2 ====")
print(table2_formatted, n = Inf)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  table2,
  file.path(TABLES_DIR, "40R_table2_absolute_risk_differences_revised_raw.csv")
)

readr::write_csv(
  table2_formatted,
  file.path(TABLES_DIR, "40R_table2_absolute_risk_differences_revised_formatted.csv")
)

saveRDS(
  table2,
  file.path(DERIVED_DIR, "40R_table2_absolute_risk_differences_revised_raw.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "40R_table2_absolute_risk_differences_revised_formatted.csv"))
message("==== REVISED TABLE 2 COMPLETE ====")