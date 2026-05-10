# =============================================================================
# 34R_fit_revised_weighted_outcome_models.R
# =============================================================================

source("R/00_setup/00N_setup_project.R")

library(dplyr)
library(broom)

message("==== FIT REVISED WEIGHTED OUTCOME MODELS ====")

# -----------------------------------------------------------------------------
# LOAD DATASETS
# -----------------------------------------------------------------------------

df_full <- readRDS("data/derived/24R_operation_dataset_full_major_noncardiac.rds")
df_abd  <- readRDS("data/derived/24R_operation_dataset_elective_abdominal_major.rds")

# -----------------------------------------------------------------------------
# MODEL FORMULA (FINAL)
# -----------------------------------------------------------------------------

form <- aki_creat_7d ~
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
# FUNCTION
# -----------------------------------------------------------------------------

fit_model <- function(df, label) {
  
  message("---- Fitting: ", label)
  
  fit <- glm(
    form,
    data = df,
    family = binomial(),
    weights = sw_1_99
  )
  
  res <- broom::tidy(fit) %>%
    mutate(
      model = label,
      OR = exp(estimate),
      CI_low = exp(estimate - 1.96 * std.error),
      CI_high = exp(estimate + 1.96 * std.error)
    )
  
  list(
    fit = fit,
    results = res
  )
}

# -----------------------------------------------------------------------------
# FIT MODELS
# -----------------------------------------------------------------------------

model_abd  <- fit_model(df_abd,  "primary_abdominal")
model_full <- fit_model(df_full, "sensitivity_full")
# Unweighted comparator for primary abdominal cohort
fit_abd_unweighted <- glm(
  form,
  data = df_abd,
  family = binomial()
)
# Unweighted comparator for sensitivity full cohort
fit_full_unweighted <- glm(
  form,
  data = df_full,
  family = binomial()
)
# -----------------------------------------------------------------------------
# COMBINE RESULTS
# -----------------------------------------------------------------------------

results <- bind_rows(
  model_abd$results,
  model_full$results
)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

saveRDS(model_abd$fit,  "data/derived/34R_model_primary_abdominal.rds")
saveRDS(model_full$fit, "data/derived/34R_model_sensitivity_full.rds")
saveRDS(
  fit_abd_unweighted,
  "data/derived/34R_model_primary_abdominal_unweighted.rds"
)
saveRDS(
  fit_full_unweighted,
  file.path(DERIVED_DIR, "34R_model_sensitivity_full_unweighted.rds")
)
write.csv(results,
          "outputs/tables/34R_model_results_combined.csv",
          row.names = FALSE)

message("==== MODELS COMPLETE ====")

# -----------------------------------------------------------------------------
# QUICK OUTPUT
# -----------------------------------------------------------------------------

results %>%
  dplyr::filter(term %in% c(
    "hypotension_time_10min",
    "map_deficit65_per_10",
    "fluid_intervals_per_5",
    "vasopressor_intervals_per_5"
  )) %>%
  dplyr::select(model, term, OR, CI_low, CI_high, p.value) %>%
  print()