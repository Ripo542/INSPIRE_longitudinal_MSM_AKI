# =============================================================================
# INSPIRE_HEMO_MSM_AKI | S1N_plot_weight_distribution.R
# Supplementary figure: stabilized weight distribution
# =============================================================================

source("R/00_setup/00N_setup_project.R")

opw <- readRDS(file.path(DERIVED_DIR, "31N_operation_final_weights.rds"))

plot_dat <- opw %>%
  dplyr::mutate(
    log10_untruncated = log10(final_sw_joint + 1e-12)
  )

p1 <- ggplot(plot_dat, aes(x = log10_untruncated)) +
  geom_histogram(bins = 80, fill = "grey55", colour = "white", linewidth = 0.15) +
  labs(
    x = expression(log[10]*"(untruncated stabilized weight)"),
    y = "Number of procedures"
  ) +
  theme_classic(base_size = 10) +
  theme(
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    axis.text = element_text(colour = "black")
  )

p2 <- ggplot(plot_dat, aes(x = final_sw_joint_trunc_1_99)) +
  geom_histogram(bins = 80, fill = "#1E40AF", colour = "white", linewidth = 0.15) +
  labs(
    x = "Stabilized weight after 1st–99th percentile truncation",
    y = "Number of procedures"
  ) +
  theme_classic(base_size = 10) +
  theme(
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
    axis.text = element_text(colour = "black")
  )

ggsave(
  file.path(FIGURES_DIR, "S1_weights_untruncated_log10.png"),
  p1, width = 5.5, height = 3.6, dpi = 600
)

ggsave(
  file.path(FIGURES_DIR, "S1_weights_truncated_1_99.png"),
  p2, width = 5.5, height = 3.6, dpi = 600
)

ggsave(
  file.path(FIGURES_DIR, "S1_weights_untruncated_log10.pdf"),
  p1, width = 5.5, height = 3.6
)

ggsave(
  file.path(FIGURES_DIR, "S1_weights_truncated_1_99.pdf"),
  p2, width = 5.5, height = 3.6
)

message("Supplementary weight figures saved.")