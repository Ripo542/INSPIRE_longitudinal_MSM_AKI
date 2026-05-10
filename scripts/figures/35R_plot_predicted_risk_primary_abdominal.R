# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 35R_plot_predicted_risk_primary_abdominal.R
# Purpose:
#   Final Figure: marginal predicted AKI risk across exposures
#   Primary cohort: elective major abdominal surgery.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(MASS)
  library(scales)
  library(readr)
})

message("==== FIGURE 35R: PRIMARY ABDOMINAL COHORT ====")

# -----------------------------------------------------------------------------
# LOAD PRIMARY DATA + MODEL
# -----------------------------------------------------------------------------

op_dat <- readRDS(file.path(DERIVED_DIR, "24R_operation_dataset_elective_abdominal_major.rds"))
fit <- readRDS(file.path(DERIVED_DIR, "34R_model_primary_abdominal.rds"))

mean_risk <- mean(op_dat$aki_creat_7d == 1, na.rm = TRUE)

message("Procedures: ", nrow(op_dat))
message("AKI risk: ", round(100 * mean_risk, 2), "%")

# -----------------------------------------------------------------------------
# EXPOSURES
# -----------------------------------------------------------------------------

exposures <- tibble::tribble(
  ~var,                          ~label,
  "hypotension_time_10min",       "Hypotension duration",
  "map_deficit65_per_10",         "MAP deficit <65 mmHg",
  "fluid_intervals_per_5",        "Fluid exposure",
  "vasopressor_intervals_per_5",  "Vasopressor exposure"
)

exposure_levels <- c(
  "Hypotension duration",
  "MAP deficit <65 mmHg",
  "Fluid exposure",
  "Vasopressor exposure"
)

# -----------------------------------------------------------------------------
# GRID: observed support, 5th–95th percentiles
# -----------------------------------------------------------------------------

ranges <- exposures %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    q05 = as.numeric(stats::quantile(op_dat[[var]], 0.05, na.rm = TRUE)),
    q95 = as.numeric(stats::quantile(op_dat[[var]], 0.95, na.rm = TRUE))
  ) %>%
  dplyr::ungroup()

grid <- ranges %>%
  dplyr::mutate(
    x = purrr::map2(q05, q95, ~ seq(.x, .y, length.out = 80))
  ) %>%
  tidyr::unnest(x)

# -----------------------------------------------------------------------------
# MARGINAL PREDICTION
# -----------------------------------------------------------------------------

predict_marginal <- function(var, value, model = fit, data = op_dat) {
  nd <- data
  nd[[var]] <- value
  mean(predict(model, newdata = nd, type = "response"), na.rm = TRUE)
}

grid <- grid %>%
  dplyr::mutate(
    y = purrr::map2_dbl(var, x, predict_marginal)
  )

# -----------------------------------------------------------------------------
# MODEL UNCERTAINTY VIA COEFFICIENT SIMULATION
# -----------------------------------------------------------------------------

set.seed(20260505)

beta_sim <- MASS::mvrnorm(
  n = 500,
  mu = coef(fit),
  Sigma = vcov(fit)
)

predict_ci <- function(var, value) {
  nd <- op_dat
  nd[[var]] <- value
  
  mm <- model.matrix(formula(fit), data = nd)
  eta <- mm %*% t(beta_sim)
  p <- plogis(eta)
  
  risk_draws <- colMeans(p, na.rm = TRUE)
  
  c(
    low = as.numeric(stats::quantile(risk_draws, 0.025, na.rm = TRUE)),
    high = as.numeric(stats::quantile(risk_draws, 0.975, na.rm = TRUE))
  )
}

ci_mat <- purrr::map2(grid$var, grid$x, predict_ci)
ci_mat <- do.call(rbind, ci_mat)

grid <- grid %>%
  dplyr::mutate(
    low = ci_mat[, "low"],
    high = ci_mat[, "high"],
    label = factor(label, levels = exposure_levels)
  )

# -----------------------------------------------------------------------------
# IQR MARKERS: P25/P75, same contrasts as Table 2
# -----------------------------------------------------------------------------

iqr_points <- exposures %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    p25 = as.numeric(stats::quantile(op_dat[[var]], 0.25, na.rm = TRUE)),
    p75 = as.numeric(stats::quantile(op_dat[[var]], 0.75, na.rm = TRUE))
  ) %>%
  dplyr::ungroup()

iqr_plot <- iqr_points %>%
  tidyr::pivot_longer(
    cols = c(p25, p75),
    names_to = "percentile",
    values_to = "x"
  ) %>%
  dplyr::mutate(
    y = purrr::map2_dbl(var, x, predict_marginal),
    percentile = dplyr::recode(
      percentile,
      p25 = "25th percentile",
      p75 = "75th percentile"
    ),
    label = factor(label, levels = exposure_levels)
  )

iqr_segments <- iqr_points %>%
  dplyr::mutate(
    y25 = purrr::map2_dbl(var, p25, predict_marginal),
    y75 = purrr::map2_dbl(var, p75, predict_marginal),
    label = factor(label, levels = exposure_levels)
  )

# -----------------------------------------------------------------------------
# EMPIRICAL SUPPORT RUG
# -----------------------------------------------------------------------------

set.seed(20260505)

rug_dat <- exposures %>%
  dplyr::mutate(
    values = purrr::map(var, ~ op_dat[[.x]])
  ) %>%
  tidyr::unnest(values) %>%
  dplyr::left_join(ranges, by = c("var", "label")) %>%
  dplyr::filter(values >= q05, values <= q95) %>%
  dplyr::group_by(var) %>%
  dplyr::mutate(.sample_prob = pmin(1, 1000 / dplyr::n())) %>%
  dplyr::ungroup() %>%
  dplyr::filter(stats::runif(dplyr::n()) < .sample_prob) %>%
  dplyr::select(-.sample_prob) %>%
  dplyr::mutate(
    label = factor(label, levels = exposure_levels)
  )
# -----------------------------------------------------------------------------
# PANEL LABELS
# -----------------------------------------------------------------------------

labeller_exposure <- function(label) {
  dplyr::case_when(
    label == "Hypotension duration" ~ "Hypotension duration\n(×10 min)",
    label == "MAP deficit <65 mmHg" ~ "MAP deficit <65 mmHg\n(per 10 units)",
    label == "Fluid exposure" ~ "Fluid exposure\n(5-min intervals)",
    label == "Vasopressor exposure" ~ "Vasopressor exposure\n(5-min intervals)",
    TRUE ~ label
  )
}

# -----------------------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------------------

p <- ggplot(grid, aes(x = x, y = y)) +
  geom_hline(
    yintercept = mean_risk,
    linetype = "dotted",
    linewidth = 0.25,
    colour = "grey65"
  ) +
  geom_vline(
    data = iqr_points,
    aes(xintercept = p25),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "grey75",
    linewidth = 0.20,
    alpha = 0.45
  ) +
  geom_vline(
    data = iqr_points,
    aes(xintercept = p75),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "grey75",
    linewidth = 0.20,
    alpha = 0.45
  ) +
  geom_ribbon(
    aes(ymin = low, ymax = high),
    fill = "#D8E2EE",
    alpha = 0.40
  ) +
  geom_line(
    colour = "#1E40AF",
    linewidth = 0.90,
    lineend = "round"
  ) +
  geom_segment(
    data = iqr_segments,
    aes(
      x = p25,
      xend = p75,
      y = y25,
      yend = y75
    ),
    inherit.aes = FALSE,
    colour = "#111827",
    linewidth = 0.38,
    alpha = 0.38,
    lineend = "round"
  ) +
  geom_point(
    data = iqr_plot,
    aes(
      x = x,
      y = y,
      shape = percentile,
      fill = percentile
    ),
    inherit.aes = FALSE,
    colour = "#111827",
    stroke = 0.60,
    size = 2.00
  ) +
  scale_shape_manual(
    name = NULL,
    values = c(
      "25th percentile" = 21,
      "75th percentile" = 24
    )
  ) +
  scale_fill_manual(
    name = NULL,
    values = c(
      "25th percentile" = "white",
      "75th percentile" = "#4B67A1"
    )
  ) +
  facet_wrap(
    ~ label,
    scales = "free_x",
    nrow = 1,
    labeller = labeller(label = labeller_exposure)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    breaks = scales::pretty_breaks(n = 4),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 4)
  ) +
  labs(
    x = NULL,
    y = "Predicted postoperative AKI risk"
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(size = 2.3, stroke = 0.65)
    ),
    fill = guide_legend(
      override.aes = list(size = 2.3, stroke = 0.65)
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "plain",
      size = 8.1,
      colour = "black",
      lineheight = 0.95
    ),
    axis.text = element_text(colour = "black", size = 8.1),
    axis.title.y = element_text(size = 9.6, margin = margin(r = 8)),
    axis.line = element_line(colour = "black", linewidth = 0.32),
    axis.ticks = element_line(colour = "black", linewidth = 0.26),
    panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.28),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1.12, "lines"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.2, colour = "black"),
    legend.key.height = unit(0.43, "cm"),
    legend.key.width = unit(0.43, "cm"),
    plot.margin = margin(5, 8, 5, 5)
  )

print(p)

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

ggsave(
  filename = file.path(FIGURES_DIR, "35R_predicted_risk_primary_abdominal.png"),
  plot = p,
  width = 8.2,
  height = 3.4,
  dpi = 600
)

ggsave(
  filename = file.path(FIGURES_DIR, "35R_predicted_risk_primary_abdominal.pdf"),
  plot = p,
  width = 8.2,
  height = 3.4
)

ggsave(
  filename = file.path(FIGURES_DIR, "35R_predicted_risk_primary_abdominal.tiff"),
  plot = p,
  width = 8.2,
  height = 3.4,
  dpi = 600,
  compression = "lzw"
)

readr::write_csv(
  grid,
  file.path(TABLES_DIR, "35R_predicted_risk_primary_abdominal.csv")
)

readr::write_csv(
  iqr_plot,
  file.path(TABLES_DIR, "35R_predicted_risk_primary_abdominal_iqr_points.csv")
)

message("Saved:")
message(file.path(FIGURES_DIR, "35R_predicted_risk_primary_abdominal.png"))
message("==== FIGURE 35R COMPLETE ====")