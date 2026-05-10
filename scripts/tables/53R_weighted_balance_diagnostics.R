# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 53R_weighted_balance_diagnostics.R
# Purpose:
#   Supplementary balance diagnostics before and after weighting.
#   Computes absolute standardized mean differences (ASMD/SMD) for:
#   1) baseline covariates
#   2) lagged time-varying history variables
#   according to current fluid and vasopressor exposure.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(stringr)
})

message("==== WEIGHTED BALANCE DIAGNOSTICS ====")

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

long <- readRDS(
  file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")
) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id)
  )

message("Rows loaded: ", nrow(long))
message("Operations loaded: ", dplyr::n_distinct(long$op_id))

# -----------------------------------------------------------------------------
# Settings
# -----------------------------------------------------------------------------

weight_col <- "final_sw_joint_trunc_1_99"

if (!weight_col %in% names(long)) {
  stop("Weight column not found in long: ", weight_col)
}

required_cols <- c(
  "revised_primary_abdominal_complete",
  "fluid_any",
  "vasopressor_any",
  weight_col
)

missing_required <- setdiff(required_cols, names(long))

if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

# -----------------------------------------------------------------------------
# Restrict to primary abdominal cohort
# -----------------------------------------------------------------------------

dat <- long %>%
  dplyr::filter(revised_primary_abdominal_complete == 1) %>%
  dplyr::mutate(
    fluid_any = as.integer(fluid_any == 1),
    vasopressor_any = as.integer(vasopressor_any == 1),
    weight_final = .data[[weight_col]]
  ) %>%
  dplyr::filter(
    !is.na(weight_final),
    is.finite(weight_final),
    weight_final > 0
  )

message("Primary abdominal rows: ", nrow(dat))
message("Primary abdominal operations: ", dplyr::n_distinct(dat$op_id))

# -----------------------------------------------------------------------------
# Construct categorical dummies for balance diagnostics
# -----------------------------------------------------------------------------

dat <- dat %>%
  dplyr::mutate(
    asa_chr = as.character(asa),
    asa_II = as.integer(asa_chr %in% c("2", "II")),
    asa_III = as.integer(asa_chr %in% c("3", "III")),
    asa_IV_V = as.integer(asa_chr %in% c("4", "5", "IV", "V")),
    
    dept_chr = as.character(department_clean),
    dept_OG = as.integer(dept_chr == "OG"),
    dept_UR = as.integer(dept_chr == "UR"),
    
    female = as.integer(female == 1),
    diabetes = as.integer(diabetes == 1),
    ckd = as.integer(ckd == 1),
    ischaemic_heart_disease = as.integer(ischaemic_heart_disease == 1),
    chronic_pulmonary_disease = as.integer(chronic_pulmonary_disease == 1),
    
    fluid_any_lag1 = as.integer(fluid_any_lag1 == 1),
    vasopressor_any_lag1 = as.integer(vasopressor_any_lag1 == 1),
    hypotension65_lag1 = as.integer(hypotension65_lag1 == 1)
  )

# -----------------------------------------------------------------------------
# Variables to assess
# -----------------------------------------------------------------------------

baseline_vars <- c(
  "age_z",
  "female",
  "asa_II",
  "asa_III",
  "asa_IV_V",
  "bmi_z",
  "baseline_creatinine_z",
  "baseline_haemoglobin",
  "baseline_albumin",
  "dept_OG",
  "dept_UR",
  "diabetes",
  "ckd",
  "ischaemic_heart_disease",
  "chronic_pulmonary_disease",
  "n_intervals"
)

timevarying_vars <- c(
  "map_lag1_z",
  "hr_lag1_z",
  "hypotension65_lag1",
  "map_deficit65_lag1_z",
  "map_change_prev_z",
  "cumulative_hypotension_intervals_prev_z",
  "cumulative_map_deficit65_prev_z",
  "cumulative_fluid_value_prev_z",
  "cumulative_vasopressor_intervals_prev_z",
  "fluid_any_lag1",
  "vasopressor_any_lag1",
  "time_elapsed_hr_z"
)

baseline_vars <- baseline_vars[baseline_vars %in% names(dat)]
timevarying_vars <- timevarying_vars[timevarying_vars %in% names(dat)]

var_labels <- c(
  age_z = "Age",
  female = "Female sex",
  asa_II = "ASA II",
  asa_III = "ASA III",
  asa_IV_V = "ASA IV–V",
  bmi_z = "Body mass index",
  baseline_creatinine_z = "Baseline creatinine",
  baseline_haemoglobin = "Baseline haemoglobin",
  baseline_albumin = "Baseline albumin",
  dept_OG = "Gynecology",
  dept_UR = "Urology",
  diabetes = "Diabetes",
  ckd = "Chronic kidney disease",
  ischaemic_heart_disease = "Ischaemic heart disease",
  chronic_pulmonary_disease = "Chronic pulmonary disease",
  n_intervals = "Number of intervals",
  
  map_lag1_z = "Prior MAP",
  hr_lag1_z = "Prior heart rate",
  hypotension65_lag1 = "Prior hypotension",
  map_deficit65_lag1_z = "Prior MAP deficit",
  map_change_prev_z = "Prior MAP change",
  cumulative_hypotension_intervals_prev_z = "Cumulative prior hypotension",
  cumulative_map_deficit65_prev_z = "Cumulative prior MAP deficit",
  cumulative_fluid_value_prev_z = "Cumulative prior fluid exposure",
  cumulative_vasopressor_intervals_prev_z = "Cumulative prior vasopressor exposure",
  fluid_any_lag1 = "Prior fluid administration",
  vasopressor_any_lag1 = "Prior vasopressor administration",
  time_elapsed_hr_z = "Elapsed time"
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  if (sum(ok) == 0) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}

weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w)
  if (sum(ok) <= 1) return(NA_real_)
  m <- weighted_mean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])
}

calc_smd_one <- function(data, treat_col, var_col, w_col = NULL) {
  x <- data[[var_col]]
  g <- data[[treat_col]]
  
  if (is.null(w_col)) {
    w <- rep(1, length(x))
  } else {
    w <- data[[w_col]]
  }
  
  ok <- !is.na(x) & !is.na(g) & !is.na(w) &
    is.finite(x) & is.finite(w) &
    g %in% c(0, 1)
  
  x <- x[ok]
  g <- g[ok]
  w <- w[ok]
  
  if (length(unique(g)) < 2) return(NA_real_)
  
  x1 <- x[g == 1]
  x0 <- x[g == 0]
  w1 <- w[g == 1]
  w0 <- w[g == 0]
  
  m1 <- weighted_mean(x1, w1)
  m0 <- weighted_mean(x0, w0)
  
  v1 <- weighted_var(x1, w1)
  v0 <- weighted_var(x0, w0)
  
  denom <- sqrt((v1 + v0) / 2)
  
  if (is.na(denom) || denom == 0) return(NA_real_)
  
  abs(m1 - m0) / denom
}

calc_balance <- function(data, treat_col, vars, panel_label, treatment_label) {
  tibble::tibble(variable = vars) %>%
    dplyr::mutate(
      treatment = treatment_label,
      panel = panel_label,
      ASMD_unweighted = purrr::map_dbl(
        variable,
        ~ calc_smd_one(data, treat_col = treat_col, var_col = .x, w_col = NULL)
      ),
      ASMD_weighted = purrr::map_dbl(
        variable,
        ~ calc_smd_one(data, treat_col = treat_col, var_col = .x, w_col = "weight_final")
      )
    )
}

# purrr check
if (!requireNamespace("purrr", quietly = TRUE)) {
  stop("Package 'purrr' is required. Install it with install.packages('purrr').")
}

# -----------------------------------------------------------------------------
# Calculate balance
# -----------------------------------------------------------------------------

balance_baseline_fluid <- calc_balance(
  dat,
  treat_col = "fluid_any",
  vars = baseline_vars,
  panel_label = "Baseline covariates",
  treatment_label = "Fluid administration"
)

balance_time_fluid <- calc_balance(
  dat,
  treat_col = "fluid_any",
  vars = timevarying_vars,
  panel_label = "Lagged time-varying history",
  treatment_label = "Fluid administration"
)

balance_baseline_vaso <- calc_balance(
  dat,
  treat_col = "vasopressor_any",
  vars = baseline_vars,
  panel_label = "Baseline covariates",
  treatment_label = "Vasopressor administration"
)

balance_time_vaso <- calc_balance(
  dat,
  treat_col = "vasopressor_any",
  vars = timevarying_vars,
  panel_label = "Lagged time-varying history",
  treatment_label = "Vasopressor administration"
)

balance_long <- dplyr::bind_rows(
  balance_baseline_fluid,
  balance_time_fluid,
  balance_baseline_vaso,
  balance_time_vaso
) %>%
  dplyr::mutate(
    label = dplyr::coalesce(unname(var_labels[variable]), variable)
  )

balance_plot <- balance_long %>%
  tidyr::pivot_longer(
    cols = c(ASMD_unweighted, ASMD_weighted),
    names_to = "model",
    values_to = "ASMD"
  ) %>%
  dplyr::mutate(
    model = dplyr::case_when(
      model == "ASMD_unweighted" ~ "Unweighted",
      model == "ASMD_weighted" ~ "Weighted",
      TRUE ~ model
    ),
    panel = factor(
      panel,
      levels = c("Baseline covariates", "Lagged time-varying history")
    ),
    treatment = factor(
      treatment,
      levels = c("Fluid administration", "Vasopressor administration")
    )
  ) %>%
  dplyr::filter(!is.na(ASMD))

# -----------------------------------------------------------------------------
# Print summary
# -----------------------------------------------------------------------------

summary_balance <- balance_plot %>%
  dplyr::group_by(treatment, panel, model) %>%
  dplyr::summarise(
    median_ASMD = median(ASMD, na.rm = TRUE),
    max_ASMD = max(ASMD, na.rm = TRUE),
    n_above_0_10 = sum(ASMD > 0.10, na.rm = TRUE),
    n_variables = dplyr::n(),
    .groups = "drop"
  )

message("==== Balance summary ====")
print(summary_balance, n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# Save tables
# -----------------------------------------------------------------------------

readr::write_csv(
  balance_long,
  file.path(TABLES_DIR, "53R_balance_diagnostics_wide.csv")
)

readr::write_csv(
  balance_plot,
  file.path(TABLES_DIR, "53R_balance_diagnostics_long.csv")
)

readr::write_csv(
  summary_balance,
  file.path(TABLES_DIR, "53R_balance_diagnostics_summary.csv")
)

saveRDS(
  list(
    balance_wide = balance_long,
    balance_long = balance_plot,
    summary = summary_balance
  ),
  file.path(DERIVED_DIR, "53R_weighted_balance_diagnostics.rds")
)

# -----------------------------------------------------------------------------
# Plot: Love plot
# -----------------------------------------------------------------------------

p <- ggplot(
  balance_plot,
  aes(
    x = ASMD,
    y = reorder(label, ASMD),
    shape = model
  )
) +
  geom_vline(
    xintercept = 0.10,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.35
  ) +
  geom_point(
    aes(fill = model),
    size = 2.2,
    colour = "black",
    stroke = 0.25,
    alpha = 0.95
  ) +
  facet_grid(
    panel ~ treatment,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_shape_manual(
    values = c(
      "Unweighted" = 21,
      "Weighted" = 24
    )
  ) +
  scale_fill_manual(
    values = c(
      "Unweighted" = "#D8E2EE",
      "Weighted" = "#1E40AF"
    )
  ) +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = c(0, 0.05, 0.10, 0.20, 0.30, 0.40),
    expand = expansion(mult = c(0.01, 0.08))
  ) +
  labs(
    x = "Absolute standardized mean difference",
    y = NULL,
    shape = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9, colour = "black"),
    axis.text = element_text(size = 7.8, colour = "black"),
    axis.title.x = element_text(size = 9.5, margin = margin(t = 7)),
    axis.line = element_line(colour = "black", linewidth = 0.30),
    axis.ticks = element_line(colour = "black", linewidth = 0.25),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.28),
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 8.5),
    panel.spacing.x = unit(1.1, "lines"),
    panel.spacing.y = unit(1.2, "lines"),
    plot.margin = margin(6, 8, 6, 6)
  )

print(p)

ggsave(
  filename = file.path(FIGURES_DIR, "S2R_weighted_balance_diagnostics_love_plot.png"),
  plot = p,
  width = 8.2,
  height = 8.6,
  dpi = 600
)

ggsave(
  filename = file.path(FIGURES_DIR, "S2R_weighted_balance_diagnostics_love_plot.pdf"),
  plot = p,
  width = 8.2,
  height = 8.6
)

ggsave(
  filename = file.path(FIGURES_DIR, "S2R_weighted_balance_diagnostics_love_plot.tiff"),
  plot = p,
  width = 8.2,
  height = 8.6,
  dpi = 600,
  compression = "lzw"
)

message("Saved:")
message(file.path(TABLES_DIR, "53R_balance_diagnostics_summary.csv"))
message(file.path(FIGURES_DIR, "S2R_weighted_balance_diagnostics_love_plot.png"))
message("==== WEIGHTED BALANCE DIAGNOSTICS COMPLETE ====")