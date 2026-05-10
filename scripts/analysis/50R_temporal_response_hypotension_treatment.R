# =============================================================================
# Temporal response: hypotension → treatment (descriptive)
# =============================================================================

source("R/00_setup/00N_setup_project.R")

library(dplyr)
library(readr)
library(tibble)

message("==== TEMPORAL RESPONSE ANALYSIS ====")

# -----------------------------------------------------------------------------
# Load longitudinal dataset
# -----------------------------------------------------------------------------

long <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds"))

# -----------------------------------------------------------------------------
# Restrict to final analytic cohorts (optional: abdominal primary)
# -----------------------------------------------------------------------------

long <- long %>%
  dplyr::filter(revised_primary_abdominal_complete == 1)

# -----------------------------------------------------------------------------
# Operation-level probabilities
# -----------------------------------------------------------------------------

op_probs <- long %>%
  dplyr::filter(!is.na(hypotension65_lag1)) %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    
    # Fluid
    p_fluid_hypo = mean(fluid_any[hypotension65_lag1 == 1], na.rm = TRUE),
    p_fluid_nohypo = mean(fluid_any[hypotension65_lag1 == 0], na.rm = TRUE),
    
    # Vasopressor
    p_vaso_hypo = mean(vasopressor_any[hypotension65_lag1 == 1], na.rm = TRUE),
    p_vaso_nohypo = mean(vasopressor_any[hypotension65_lag1 == 0], na.rm = TRUE),
    
    # Combined
    p_both_hypo = mean(
      (fluid_any == 1 & vasopressor_any == 1)[hypotension65_lag1 == 1],
      na.rm = TRUE
    ),
    p_both_nohypo = mean(
      (fluid_any == 1 & vasopressor_any == 1)[hypotension65_lag1 == 0],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

summary_tbl <- tibble::tibble(
  
  Metric = c(
    "Fluid administration",
    "Vasopressor administration",
    "Combined fluid + vasopressor"
  ),
  
  Hypotension_prev = c(
    mean(op_probs$p_fluid_hypo, na.rm = TRUE),
    mean(op_probs$p_vaso_hypo, na.rm = TRUE),
    mean(op_probs$p_both_hypo, na.rm = TRUE)
  ),
  
  No_hypotension_prev = c(
    mean(op_probs$p_fluid_nohypo, na.rm = TRUE),
    mean(op_probs$p_vaso_nohypo, na.rm = TRUE),
    mean(op_probs$p_both_nohypo, na.rm = TRUE)
  )
) %>%
  dplyr::mutate(
    Absolute_difference = Hypotension_prev - No_hypotension_prev
  )

print(summary_tbl)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  summary_tbl,
  file.path(TABLES_DIR, "50R_temporal_response_summary.csv")
)

saveRDS(
  op_probs,
  file.path(DERIVED_DIR, "50R_temporal_response_op_level.rds")
)

message("Saved.")