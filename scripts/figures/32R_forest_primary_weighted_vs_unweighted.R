# =============================================================================
# 32R_forest_primary_weighted_vs_unweighted.R
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(broom)
  library(readr)
})

fit_w <- readRDS(file.path(DERIVED_DIR, "34R_model_primary_abdominal.rds"))
fit_u <- readRDS(file.path(DERIVED_DIR, "34R_model_primary_abdominal_unweighted.rds"))

terms_keep <- c(
  "hypotension_time_10min",
  "map_deficit65_per_10",
  "fluid_intervals_per_5",
  "vasopressor_intervals_per_5"
)

term_labels <- c(
  hypotension_time_10min = "Hypotension duration\n(per 10 min)",
  map_deficit65_per_10 = "MAP deficit <65 mmHg\n(per 10 units)",
  fluid_intervals_per_5 = "Fluid exposure\n(5-min intervals)",
  vasopressor_intervals_per_5 = "Vasopressor exposure\n(5-min intervals)"
)

plot_dat <- dplyr::bind_rows(
  broom::tidy(fit_w) %>%
    dplyr::filter(term %in% terms_keep) %>%
    dplyr::mutate(model = "Weighted longitudinal model"),
  broom::tidy(fit_u) %>%
    dplyr::filter(term %in% terms_keep) %>%
    dplyr::mutate(model = "Unweighted comparator")
) %>%
  dplyr::mutate(
    or = exp(estimate),
    low = exp(estimate - 1.96 * std.error),
    high = exp(estimate + 1.96 * std.error),
    exposure = dplyr::recode(term, !!!term_labels),
    exposure = factor(
      exposure,
      levels = rev(c(
        "Hypotension duration\n(per 10 min)",
        "MAP deficit <65 mmHg\n(per 10 units)",
        "Fluid exposure\n(5-min intervals)",
        "Vasopressor exposure\n(5-min intervals)"
      ))
    ),
    model = factor(
      model,
      levels = c("Weighted longitudinal model", "Unweighted comparator")
    )
  )

p <- ggplot(
  plot_dat,
  aes(
    x = or,
    y = exposure,
    colour = model,
    group = model
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey55",
    linewidth = 0.35
  ) +
  geom_errorbarh(
    aes(xmin = low, xmax = high),
    height = 0.12,
    linewidth = 0.55,
    position = position_dodge(width = 0.38)
  ) +
  geom_point(
    size = 2.4,
    position = position_dodge(width = 0.38)
  ) +
  scale_colour_manual(
    values = c(
      "Weighted longitudinal model" = "#1E40AF",
      "Unweighted comparator" = "#B8C0CC"
    )
  ) +
  scale_x_log10(
    breaks = c(0.90, 0.95, 1.00, 1.05, 1.10, 1.20),
    limits = c(0.90, 1.20)
  ) +
  labs(
    x = "Odds ratio for postoperative AKI",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8.7, colour = "black"),
    axis.text.x = element_text(size = 8.5, colour = "black"),
    axis.title.x = element_text(size = 9.5, margin = margin(t = 6)),
    axis.line = element_line(colour = "black", linewidth = 0.32),
    axis.ticks = element_line(colour = "black", linewidth = 0.28),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.30),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8.5),
    legend.key.width = unit(0.7, "cm"),
    plot.margin = margin(5, 8, 5, 5)
  )

print(p)

ggsave(
  file.path(FIGURES_DIR, "32R_forest_primary_weighted_vs_unweighted.png"),
  p,
  width = 7.0,
  height = 4.0,
  dpi = 600
)

ggsave(
  file.path(FIGURES_DIR, "32R_forest_primary_weighted_vs_unweighted.pdf"),
  p,
  width = 7.0,
  height = 4.0
)

ggsave(
  file.path(FIGURES_DIR, "32R_forest_primary_weighted_vs_unweighted.tiff"),
  p,
  width = 7.0,
  height = 4.0,
  dpi = 600,
  compression = "lzw"
)

readr::write_csv(
  plot_dat,
  file.path(TABLES_DIR, "32R_forest_primary_weighted_vs_unweighted.csv")
)