# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 40R_build_table2_weighted_vs_unweighted_ard.R
# Purpose:
#   Build Table 2 for the primary elective abdominal cohort:
#   weighted longitudinal model vs unweighted comparator.
#   Absolute risk differences estimated across P25–P75 exposure contrasts.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== TABLE 2: WEIGHTED VS UNWEIGHTED ARD ====")

# -----------------------------------------------------------------------------
# Load primary abdominal dataset and models
# -----------------------------------------------------------------------------

dat <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_full_major_noncardiac.rds"))

fit_weighted <- readRDS(file.path(DERIVED_DIR, "34R_model_sensitivity_full.rds"))
fit_unweighted <- readRDS(file.path(DERIVED_DIR, "34R_model_sensitivity_full_unweighted.rds"))

dat <- dat %>%
  dplyr::mutate(
    asa_final = as.factor(asa_final),
    department_final = as.factor(department_final),
    diabetes = as.integer(diabetes),
    ckd = as.integer(ckd),
    ischaemic_heart_disease = as.integer(ischaemic_heart_disease),
    chronic_pulmonary_disease = as.integer(chronic_pulmonary_disease)
  )

message("Procedures: ", nrow(dat))
message("AKI events: ", sum(dat$aki_creat_7d == 1, na.rm = TRUE))
message("AKI risk: ", round(100 * mean(dat$aki_creat_7d == 1, na.rm = TRUE), 2), "%")

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
# ARD function
# -----------------------------------------------------------------------------

compute_ard_boot <- function(dat, fit, model_label, var, exposure_label,
                             weighted = TRUE, n_boot = 1000, seed = 20260505) {
  
  message("Computing: ", model_label, " | ", exposure_label)
  
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
    boot_dat$department_final <- factor(
      boot_dat$department_final,
      levels = levels(dat$department_final)
    )
    
    fit_boot <- tryCatch(
      {
        if (weighted) {
          glm(
            model_formula,
            data = boot_dat,
            weights = sw_1_99,
            family = quasibinomial()
          )
        } else {
          glm(
            model_formula,
            data = boot_dat,
            family = quasibinomial()
          )
        }
      },
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
    Model = model_label,
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
# Build long table
# -----------------------------------------------------------------------------

exposures <- list(
  list(var = "hypotension_time_10min", exposure = "Hypotension duration"),
  list(var = "map_deficit65_per_10", exposure = "MAP deficit <65 mmHg"),
  list(var = "fluid_intervals_per_5", exposure = "Fluid exposure"),
  list(var = "vasopressor_intervals_per_5", exposure = "Vasopressor exposure")
)

table2_long <- dplyr::bind_rows(
  lapply(seq_along(exposures), function(i) {
    compute_ard_boot(
      dat = dat,
      fit = fit_weighted,
      model_label = "Weighted longitudinal model",
      var = exposures[[i]]$var,
      exposure_label = exposures[[i]]$exposure,
      weighted = TRUE,
      n_boot = 1000,
      seed = 20260500 + i
    )
  }),
  lapply(seq_along(exposures), function(i) {
    compute_ard_boot(
      dat = dat,
      fit = fit_unweighted,
      model_label = "Unweighted comparator",
      var = exposures[[i]]$var,
      exposure_label = exposures[[i]]$exposure,
      weighted = FALSE,
      n_boot = 1000,
      seed = 20260600 + i
    )
  })
)

# -----------------------------------------------------------------------------
# Format long version
# -----------------------------------------------------------------------------

table2_long_formatted <- table2_long %>%
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
    Model,
    Exposure,
    Contrast,
    Low_value,
    High_value,
    Predicted_risk_low,
    Predicted_risk_high,
    Absolute_risk_difference,
    Bootstrap_iterations
  )

# -----------------------------------------------------------------------------
# Wide version for manuscript
# -----------------------------------------------------------------------------

table2_wide <- table2_long %>%
  dplyr::mutate(
    ARD_label = sprintf(
      "%+.2f%% (%.2f to %.2f)",
      100 * Absolute_risk_difference,
      100 * CI_low,
      100 * CI_high
    ),
    Risk_low_label = sprintf("%.1f%%", 100 * Predicted_risk_low),
    Risk_high_label = sprintf("%.1f%%", 100 * Predicted_risk_high),
    Low_value = round(Low_value, 2),
    High_value = round(High_value, 2)
  ) %>%
  dplyr::select(
    Exposure,
    Contrast,
    Low_value,
    High_value,
    Model,
    Risk_low_label,
    Risk_high_label,
    ARD_label
  ) %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = c(Risk_low_label, Risk_high_label, ARD_label),
    names_glue = "{Model}_{.value}"
  ) %>%
  dplyr::rename(
    Weighted_predicted_risk_low = `Weighted longitudinal model_Risk_low_label`,
    Weighted_predicted_risk_high = `Weighted longitudinal model_Risk_high_label`,
    Weighted_ARD = `Weighted longitudinal model_ARD_label`,
    Unweighted_predicted_risk_low = `Unweighted comparator_Risk_low_label`,
    Unweighted_predicted_risk_high = `Unweighted comparator_Risk_high_label`,
    Unweighted_ARD = `Unweighted comparator_ARD_label`
  )

message("==== Table 2 long ====")
print(table2_long_formatted, n = Inf)

message("==== Table 2 wide ====")
print(table2_wide, n = Inf)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  table2_long,
  file.path(TABLES_DIR, "40S_table2_full_weighted_vs_unweighted_ard_raw_long.csv")
)

readr::write_csv(
  table2_long_formatted,
  file.path(TABLES_DIR, "40S_table2_full_weighted_vs_unweighted_ard_formatted_long.csv")
)

readr::write_csv(
  table2_wide,
  file.path(TABLES_DIR, "40S_table2_full_weighted_vs_unweighted_ard_formatted_wide.csv")
)

saveRDS(
  table2_long,
  file.path(DERIVED_DIR, "40S_table2_full_weighted_vs_unweighted_ard_raw_long.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "40S_table2_full_weighted_vs_unweighted_ard_formatted_wide.csv"))
message("==== TABLE S (FULL COHORT) COMPLETE ====")