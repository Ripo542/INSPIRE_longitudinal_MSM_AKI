# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 34N_plot_final_forest.R
# Final publication-ready forest plot with visual shift between models
# =============================================================================

source("R/00_setup/00N_setup_project.R")

message("==== PLOT FINAL FOREST FIGURE ====")

res <- readr::read_csv(
  file.path(TABLES_DIR, "34N_final_outcome_model_main_terms.csv"),
  show_col_types = FALSE
)

plot_dat <- res %>%
  dplyr::filter(
    model %in% c("primary_weighted_1_99", "unweighted_comparator"),
    term %in% c(
      "hypotension_time_10min",
      "map_deficit65_per_10",
      "fluid_intervals_per_5",
      "vasopressor_intervals_per_5"
    )
  ) %>%
  dplyr::mutate(
    exposure = dplyr::case_when(
      term == "hypotension_time_10min" ~ "Hypotension duration\n(per 10 min)",
      term == "map_deficit65_per_10" ~ "MAP deficit <65 mmHg\n(per 10 units)",
      term == "fluid_intervals_per_5" ~ "Fluid exposure\n(per 5-min interval)",
      term == "vasopressor_intervals_per_5" ~ "Vasopressor exposure\n(per 5-min interval)",
      TRUE ~ term
    ),
    model_label = dplyr::case_when(
      model == "primary_weighted_1_99" ~ "Weighted longitudinal model",
      model == "unweighted_comparator" ~ "Unweighted comparator",
      TRUE ~ model
    ),
    exposure = factor(
      exposure,
      levels = rev(c(
        "Hypotension duration\n(per 10 min)",
        "MAP deficit <65 mmHg\n(per 10 units)",
        "Fluid exposure\n(per 5-min interval)",
        "Vasopressor exposure\n(per 5-min interval)"
      ))
    ),
    model_label = factor(
      model_label,
      levels = c("Weighted longitudinal model", "Unweighted comparator")
    )
  )

# ---- Data for visual shift lines: unweighted -> weighted ----
shift_dat <- plot_dat %>%
  dplyr::select(term, exposure, model_label, or) %>%
  tidyr::pivot_wider(
    names_from = model_label,
    values_from = or
  ) %>%
  dplyr::mutate(
    x_start = `Unweighted comparator`,
    x_end = `Weighted longitudinal model`
  )

# ---- Plot ----
p <- ggplot() +
  geom_vline(
    xintercept = 1,
    linewidth = 0.30,
    linetype = "dashed",
    colour = "grey55"
  ) +
  geom_segment(
    data = shift_dat,
    aes(
      x = x_start,
      xend = x_end,
      y = exposure,
      yend = exposure
    ),
    colour = "grey85",
    linewidth = 0.60,
    alpha = 0.60,
    lineend = "round"
  ) +
  geom_errorbarh(
    data = plot_dat,
    aes(
      x = or,
      y = exposure,
      xmin = ci_low,
      xmax = ci_high,
      colour = model_label,
      linewidth = model_label
    ),
    height = 0.13,
    position = position_dodge(width = 0.46)
  ) +
  geom_point(
    data = plot_dat,
    aes(
      x = or,
      y = exposure,
      colour = model_label,
      size = model_label
    ),
    position = position_dodge(width = 0.46)
  ) +
  scale_colour_manual(
    values = c(
      "Weighted longitudinal model" = "#1E40AF",
      "Unweighted comparator" = "#B0B7C3"
    )
  ) +
  scale_size_manual(
    values = c(
      "Weighted longitudinal model" = 2.9,
      "Unweighted comparator" = 2.2
    ),
    guide = "none"
  ) +
  scale_linewidth_manual(
    values = c(
      "Weighted longitudinal model" = 0.65,
      "Unweighted comparator" = 0.45
    ),
    guide = "none"
  ) +
  scale_x_log10(
    breaks = c(0.95, 1.00, 1.05, 1.10, 1.20),
    limits = c(0.94, 1.22)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Odds ratio for postoperative AKI",
    y = NULL,
    colour = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9.5, colour = "black"),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.title.x = element_text(size = 10, margin = margin(t = 8)),
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    legend.key.width = unit(1.0, "cm"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.30),
    plot.margin = margin(5, 10, 5, 5)
  )

print(p)

# ---- Save ----
ggsave(
  filename = file.path(FIGURES_DIR, "34N_final_forest_shift_pub.png"),
  plot = p,
  width = 5.8,
  height = 3.8,
  dpi = 600
)

ggsave(
  filename = file.path(FIGURES_DIR, "34N_final_forest_shift_pub.pdf"),
  plot = p,
  width = 5.8,
  height = 3.8
)

ggsave(
  filename = file.path(FIGURES_DIR, "34N_final_forest_shift_pub.tiff"),
  plot = p,
  width = 5.8,
  height = 3.8,
  dpi = 600,
  compression = "lzw"
)

message("Saved:")
message(file.path(FIGURES_DIR, "34N_final_forest_shift_pub.png"))
message(file.path(FIGURES_DIR, "34N_final_forest_shift_pub.pdf"))
message(file.path(FIGURES_DIR, "34N_final_forest_shift_pub.tiff"))