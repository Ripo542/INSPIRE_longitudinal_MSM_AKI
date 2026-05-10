# =============================================================================
# INSPIRE_HEMO_MSM_AKI | F4R_treatment_probability_by_prior_map.R
# Purpose:
#   Figure: probability of treatment in current interval according to prior MAP.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(mgcv)
  library(scales)
})

message("==== FIGURE: TREATMENT PROBABILITY BY PRIOR MAP ====")

# -----------------------------------------------------------------------------
# Load longitudinal dataset
# -----------------------------------------------------------------------------

long <- readRDS(
  file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")
)

# -----------------------------------------------------------------------------
# Restrict to primary abdominal cohort
# -----------------------------------------------------------------------------

dat <- long %>%
  dplyr::filter(revised_primary_abdominal_complete == 1) %>%
  dplyr::filter(
    !is.na(map_lag1),
    is.finite(map_lag1),
    map_lag1 >= 45,
    map_lag1 <= 110
  ) %>%
  dplyr::mutate(
    fluid_any = as.integer(fluid_any == 1),
    vasopressor_any = as.integer(vasopressor_any == 1),
    both_any = as.integer(fluid_any == 1 & vasopressor_any == 1)
  )

message("Rows used: ", nrow(dat))
message("Operations used: ", dplyr::n_distinct(dat$op_id))

# -----------------------------------------------------------------------------
# Fit GAMs
# -----------------------------------------------------------------------------

fit_fluid <- mgcv::gam(
  fluid_any ~ s(map_lag1, k = 8),
  data = dat,
  family = binomial(link = "logit"),
  method = "REML"
)

fit_vaso <- mgcv::gam(
  vasopressor_any ~ s(map_lag1, k = 8),
  data = dat,
  family = binomial(link = "logit"),
  method = "REML"
)

fit_both <- mgcv::gam(
  both_any ~ s(map_lag1, k = 8),
  data = dat,
  family = binomial(link = "logit"),
  method = "REML"
)

# -----------------------------------------------------------------------------
# Prediction grid
# -----------------------------------------------------------------------------

grid <- tibble::tibble(
  map_lag1 = seq(
    quantile(dat$map_lag1, 0.01, na.rm = TRUE),
    quantile(dat$map_lag1, 0.99, na.rm = TRUE),
    length.out = 250
  )
)

predict_prob <- function(fit, grid, label) {
  pred <- predict(
    fit,
    newdata = grid,
    type = "link",
    se.fit = TRUE
  )
  
  tibble::tibble(
    map_lag1 = grid$map_lag1,
    outcome = label,
    fit = plogis(pred$fit),
    low = plogis(pred$fit - 1.96 * pred$se.fit),
    high = plogis(pred$fit + 1.96 * pred$se.fit)
  )
}

pred_dat <- dplyr::bind_rows(
  predict_prob(fit_fluid, grid, "Fluid administration"),
  predict_prob(fit_vaso, grid, "Vasopressor administration"),
  predict_prob(fit_both, grid, "Fluids + vasopressors")
) %>%
  dplyr::mutate(
    outcome = factor(
      outcome,
      levels = c(
        "Fluid administration",
        "Vasopressor administration",
        "Fluids + vasopressors"
      )
    )
  )

# -----------------------------------------------------------------------------
# Rug data: sampled for plotting only
# -----------------------------------------------------------------------------

set.seed(542)

n_rug <- min(12000, nrow(dat))

rug_dat <- dat %>%
  dplyr::select(map_lag1) %>%
  dplyr::slice_sample(n = n_rug)
# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

p <- ggplot(pred_dat, aes(x = map_lag1, y = fit)) +
  geom_ribbon(
    aes(ymin = low, ymax = high),
    fill = "#D8E2EE",
    alpha = 0.55
  ) +
  geom_line(
    colour = "#1E40AF",
    linewidth = 0.85,
    lineend = "round"
  ) +
  geom_vline(
    xintercept = 65,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.32
  ) +
  geom_rug(
    data = rug_dat,
    aes(x = map_lag1),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.018,
    colour = "grey45",
    linewidth = 0.15
  ) +
  facet_wrap(
    ~ outcome,
    ncol = 1,
    scales = "free_y",
    strip.position = "top"
  ) +
  scale_x_continuous(
    breaks = seq(50, 110, by = 10),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = "MAP in the preceding 5-min interval (mmHg)",
    y = "Probability of treatment in the current interval"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 9,
      face = "plain",
      colour = "black",
      margin = margin(b = 4)
    ),
    axis.text = element_text(size = 8.5, colour = "black"),
    axis.title.x = element_text(size = 9.5, margin = margin(t = 7)),
    axis.title.y = element_text(size = 9.5, margin = margin(r = 7)),
    axis.line = element_line(colour = "black", linewidth = 0.32),
    axis.ticks = element_line(colour = "black", linewidth = 0.28),
    panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.28),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1.05, "lines"),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 5)
  )

print(p)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

ggsave(
  filename = file.path(FIGURES_DIR, "F4R_treatment_probability_by_prior_map.png"),
  plot = p,
  width = 6.3,
  height = 6.2,
  dpi = 600
)

ggsave(
  filename = file.path(FIGURES_DIR, "F4R_treatment_probability_by_prior_map.pdf"),
  plot = p,
  width = 6.3,
  height = 6.2
)

ggsave(
  filename = file.path(FIGURES_DIR, "F4R_treatment_probability_by_prior_map.tiff"),
  plot = p,
  width = 6.3,
  height = 6.2,
  dpi = 600,
  compression = "lzw"
)

readr::write_csv(
  pred_dat,
  file.path(TABLES_DIR, "F4R_treatment_probability_by_prior_map_predictions.csv")
)

message("Saved:")
message(file.path(FIGURES_DIR, "F4R_treatment_probability_by_prior_map.png"))
message("==== FIGURE COMPLETE ====")