# =============================================================================
# 52R_vasopressor_type_descriptive.R
# Purpose:
#   Descriptive analysis of vasopressor exposure by type (ephedrine,
#   phenylephrine, norepinephrine) within the primary abdominal cohort.
#   No modelling. No assumptions on coding beyond existing flags.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

message("==== VASOPRESSOR TYPE DESCRIPTIVE ====")

# -----------------------------------------------------------------------------
# Load longitudinal dataset
# -----------------------------------------------------------------------------

long <- readRDS(
  file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")
) %>%
  mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id)
  )

message("Rows: ", nrow(long))
message("Operations: ", n_distinct(long$op_id))

# -----------------------------------------------------------------------------
# Check available vasopressor columns
# -----------------------------------------------------------------------------

vaso_cols <- grep("eph|phe|nepi|vaso", names(long), value = TRUE, ignore.case = TRUE)
message("Detected vaso-related columns:")
print(vaso_cols)

# -----------------------------------------------------------------------------
# Safety: require basic variables
# -----------------------------------------------------------------------------

required_cols <- c("op_id", "interval_id", "revised_primary_abdominal_complete")
missing_req <- setdiff(required_cols, names(long))

if (length(missing_req) > 0) {
  stop("Missing required columns: ", paste(missing_req, collapse = ", "))
}

# -----------------------------------------------------------------------------
# Restrict to primary abdominal cohort
# -----------------------------------------------------------------------------

df <- long %>%
  filter(revised_primary_abdominal_complete == 1)

message("Primary cohort operations: ", n_distinct(df$op_id))

# -----------------------------------------------------------------------------
# Detect vasopressor variables (flexible, no assumptions)
# -----------------------------------------------------------------------------

detect_var <- function(pattern) {
  vars <- grep(pattern, names(df), value = TRUE, ignore.case = TRUE)
  if (length(vars) == 0) return(NULL)
  vars[1]
}

eph_var  <- detect_var("eph")
phe_var  <- detect_var("phe")
nepi_var <- detect_var("nepi")

message("Detected variables:")
message("Ephedrine: ", eph_var)
message("Phenylephrine: ", phe_var)
message("Norepinephrine: ", nepi_var)

# -----------------------------------------------------------------------------
# Build interval-level exposure flags
# -----------------------------------------------------------------------------

df <- df %>%
  mutate(
    eph_any  = if (!is.null(eph_var))  as.integer(.data[[eph_var]] > 0) else 0L,
    phe_any  = if (!is.null(phe_var))  as.integer(.data[[phe_var]] > 0) else 0L,
    nepi_any = if (!is.null(nepi_var)) as.integer(.data[[nepi_var]] > 0) else 0L
  ) %>%
  mutate(
    vaso_any = as.integer(eph_any == 1 | phe_any == 1 | nepi_any == 1)
  )

# -----------------------------------------------------------------------------
# Operation-level summary (ANY exposure)
# -----------------------------------------------------------------------------

op_level <- df %>%
  group_by(op_id) %>%
  summarise(
    eph_any_op  = max(eph_any, na.rm = TRUE),
    phe_any_op  = max(phe_any, na.rm = TRUE),
    nepi_any_op = max(nepi_any, na.rm = TRUE),
    vaso_any_op = max(vaso_any, na.rm = TRUE),
    .groups = "drop"
  )

op_summary <- tibble(
  Metric = c(
    "Any vasopressor",
    "Ephedrine",
    "Phenylephrine",
    "Norepinephrine"
  ),
  Procedures = c(
    sum(op_level$vaso_any_op),
    sum(op_level$eph_any_op),
    sum(op_level$phe_any_op),
    sum(op_level$nepi_any_op)
  ),
  Percent = 100 * Procedures / n_distinct(op_level$op_id)
)

# -----------------------------------------------------------------------------
# Interval-level exposure (density)
# -----------------------------------------------------------------------------

interval_summary <- tibble(
  Metric = c(
    "Any vasopressor",
    "Ephedrine",
    "Phenylephrine",
    "Norepinephrine"
  ),
  Intervals = c(
    sum(df$vaso_any),
    sum(df$eph_any),
    sum(df$phe_any),
    sum(df$nepi_any)
  ),
  Percent = 100 * Intervals / nrow(df)
)

# -----------------------------------------------------------------------------
# Co-exposure (important clinically)
# -----------------------------------------------------------------------------

combo_summary <- df %>%
  mutate(
    combo = case_when(
      eph_any + phe_any + nepi_any == 0 ~ "None",
      eph_any + phe_any + nepi_any == 1 ~ "Single agent",
      TRUE ~ "Multiple agents"
    )
  ) %>%
  count(combo) %>%
  mutate(
    Percent = 100 * n / sum(n)
  )

# -----------------------------------------------------------------------------
# Print
# -----------------------------------------------------------------------------

message("==== Operation-level exposure ====")
print(op_summary, n = Inf)

message("==== Interval-level exposure ====")
print(interval_summary, n = Inf)

message("==== Co-exposure ====")
print(combo_summary, n = Inf)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

write_csv(op_summary,
          file.path(TABLES_DIR, "52R_vasopressor_type_op_level.csv")
)

write_csv(interval_summary,
          file.path(TABLES_DIR, "52R_vasopressor_type_interval_level.csv")
)

write_csv(combo_summary,
          file.path(TABLES_DIR, "52R_vasopressor_type_combo.csv")
)

message("Saved tables to TABLES_DIR")
message("==== DONE ====")