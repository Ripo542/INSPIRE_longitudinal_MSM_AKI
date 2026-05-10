# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 40N_build_table2_absolute_risk_differences.R
# Purpose:
#   Build Table 2: absolute risk differences using cluster bootstrap with refit.
#   Unit of resampling: operation.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== BUILD TABLE 2: ABSOLUTE RISK DIFFERENCES WITH BOOTSTRAP REFIT ====")

# -----------------------------------------------------------------------------
# Load final operation-level dataset and primary model
# -----------------------------------------------------------------------------

op_dat <- readRDS(file.path(DERIVED_DIR, "34N_final_operation_level_dataset.rds"))
fit_primary <- readRDS(file.path(DERIVED_DIR, "34N_primary_weighted_model.rds"))

required_vars <- c(
  "aki_creat_7d",
  "hypotension_time_10min",
  "map_deficit65_per_10",
  "fluid_intervals_per_5",
  "vasopressor_intervals_per_5",
  "age_z",
  "female",
  "asa_final",
  "bmi_z",
  "baseline_creatinine_z",
  "department_final",
  "n_intervals_z",
  "sw_1_99"
)

missing_vars <- setdiff(required_vars, names(op_dat))
if (length(missing_vars) > 0) {
  stop(
    "Missing required variables in 34N_final_operation_level_dataset.rds: ",
    paste(missing_vars, collapse = ", ")
  )
}

op_dat <- op_dat %>%
  dplyr::mutate(
    asa_final = as.factor(asa_final),
    department_final = as.factor(department_final)
  )

message("Operations: ", nrow(op_dat))
message("AKI events: ", sum(op_dat$aki_creat_7d == 1, na.rm = TRUE))
message("Observed AKI risk: ", round(100 * mean(op_dat$aki_creat_7d == 1, na.rm = TRUE), 2), "%")

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
  department_final +
  n_intervals_z

# -----------------------------------------------------------------------------
# Function: ARD with bootstrap refit
# -----------------------------------------------------------------------------

compute_ard_boot <- function(var, label, n_boot = 200, seed = 20260505) {
  
  message("Computing ARD for: ", label)
  
  p25 <- as.numeric(stats::quantile(op_dat[[var]], 0.25, na.rm = TRUE))
  p75 <- as.numeric(stats::quantile(op_dat[[var]], 0.75, na.rm = TRUE))
  
  low_dat <- op_dat
  high_dat <- op_dat
  
  low_dat[[var]] <- p25
  high_dat[[var]] <- p75
  
  p_low <- predict(fit_primary, newdata = low_dat, type = "response")
  p_high <- predict(fit_primary, newdata = high_dat, type = "response")
  
  ard <- mean(p_high - p_low, na.rm = TRUE)
  
  set.seed(seed)
  
  boot <- numeric(n_boot)
  boot[] <- NA_real_
  
  for (b in seq_len(n_boot)) {
    
    idx <- sample(seq_len(nrow(op_dat)), replace = TRUE)
    boot_dat <- op_dat[idx, , drop = FALSE]
    
    # Preserve factor levels to avoid prediction/model matrix failures
    boot_dat$asa_final <- factor(boot_dat$asa_final, levels = levels(op_dat$asa_final))
    boot_dat$department_final <- factor(
      boot_dat$department_final,
      levels = levels(op_dat$department_final)
    )
    
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
  
  if (length(boot) < n_boot * 0.8) {
    warning(
      "Fewer than 80% bootstrap iterations succeeded for ",
      label,
      ": ",
      length(boot),
      " / ",
      n_boot
    )
  }
  
  tibble::tibble(
    Exposure = label,
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

## -----------------------------------------------------------------------------
# Build table
# -----------------------------------------------------------------------------

table2 <- dplyr::bind_rows(
  compute_ard_boot(
    "map_deficit65_per_10",
    "MAP deficit <65 mmHg",
    n_boot = 1000,
    seed = 20260505
  ),
  compute_ard_boot(
    "fluid_intervals_per_5",
    "Fluid exposure",
    n_boot = 1000,
    seed = 20260506
  ),
  compute_ard_boot(
    "vasopressor_intervals_per_5",
    "Vasopressor exposure",
    n_boot = 1000,
    seed = 20260507
  )
)

table2_formatted <- table2 %>%
  dplyr::mutate(
    Low_value = round(Low_value, 2),
    High_value = round(High_value, 2),
    Predicted_risk_low = sprintf("%.1f%%", 100 * Predicted_risk_low),
    Predicted_risk_high = sprintf("%.1f%%", 100 * Predicted_risk_high),
    Absolute_risk_difference = sprintf(
      "%.2f%% (%.2f to %.2f)",
      100 * Absolute_risk_difference,
      100 * CI_low,
      100 * CI_high
    )
  ) %>%
  dplyr::select(
    Exposure,
    Contrast,
    Low_value,
    High_value,
    Predicted_risk_low,
    Predicted_risk_high,
    Absolute_risk_difference,
    Bootstrap_iterations
  )

message("==== Table 2 formatted ====")
print(table2_formatted)
# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  table2,
  file.path(TABLES_DIR, "40N_table2_absolute_risk_differences_bootstrap_raw.csv")
)

readr::write_csv(
  table2_formatted,
  file.path(TABLES_DIR, "40N_table2_absolute_risk_differences_bootstrap_formatted.csv")
)

saveRDS(
  table2,
  file.path(DERIVED_DIR, "40N_table2_absolute_risk_differences_bootstrap_raw.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "40N_table2_absolute_risk_differences_bootstrap_formatted.csv"))
message("==== TABLE 2 COMPLETE ====")