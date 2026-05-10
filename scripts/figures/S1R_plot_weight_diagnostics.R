# =============================================================================
# INSPIRE_HEMO_MSM_AKI | S1R_plot_weight_diagnostics.R
# Purpose:
#   Supplementary Figure S1: distribution of final stabilised weights
#   after 1st–99th percentile truncation.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
})

message("==== SUPPLEMENTARY FIGURE S1: WEIGHT DISTRIBUTION ====")

long <- readRDS(file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds"))

extract_weights <- function(data, flag, label) {
  data %>%
    dplyr::filter(.data[[flag]] == 1) %>%
    dplyr::arrange(op_id, interval_id) %>%
    dplyr::group_by(op_id) %>%
    dplyr::summarise(
      weight = dplyr::last(final_sw_joint_trunc_1_99),
      .groups = "drop"
    ) %>%
    dplyr::mutate(cohort = label)
}

weights_plot <- dplyr::bind_rows(
  extract_weights(
    long,
    "revised_primary_abdominal_complete",
    "Primary abdominal cohort"
  ),
  extract_weights(
    long,
    "revised_primary_full_complete",
    "Full major non-cardiac cohort"
  )
) %>%
  dplyr::filter(is.finite(weight), weight > 0)

summary_lines <- weights_plot %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(
    median = stats::median(weight, na.rm = TRUE),
    p99 = as.numeric(stats::quantile(weight, 0.99, na.rm = TRUE)),
    .groups = "drop"
  )

message("==== Weight summary used in plot ====")
print(summary_lines)

p <- ggplot(weights_plot, aes(x = weight)) +
  geom_histogram(
    bins = 70,
    fill = "#D8E2EE",
    colour = "white",
    linewidth = 0.15
  ) +
  geom_vline(
    data = summary_lines,
    aes(xintercept = median),
    linetype = "dashed",
    colour = "#1E40AF",
    linewidth = 0.45
  ) +
  geom_vline(
    data = summary_lines,
    aes(xintercept = p99),
    linetype = "dotted",
    colour = "#111827",
    linewidth = 0.45
  ) +
  facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
  scale_x_log10(
    breaks = c(0.01, 0.1, 0.5, 1, 3, 10, 30),
    labels = c("0.01", "0.1", "0.5", "1", "3", "10", "30")
  ) +
  labs(
    x = "Final stabilised inverse probability weight (log scale)",
    y = "Number of procedures"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "plain", size = 9, colour = "black"),
    axis.text = element_text(size = 8.5, colour = "black"),
    axis.title = element_text(size = 9.5, colour = "black"),
    axis.line = element_line(colour = "black", linewidth = 0.32),
    axis.ticks = element_line(colour = "black", linewidth = 0.28),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(1.0, "lines"),
    plot.margin = margin(5, 8, 5, 5)
  )

print(p)

ggsave(
  filename = file.path(FIGURES_DIR, "S1R_weight_distribution_trunc_1_99.png"),
  plot = p,
  width = 6.2,
  height = 4.4,
  dpi = 600
)

ggsave(
  filename = file.path(FIGURES_DIR, "S1R_weight_distribution_trunc_1_99.pdf"),
  plot = p,
  width = 6.2,
  height = 4.4
)

ggsave(
  filename = file.path(FIGURES_DIR, "S1R_weight_distribution_trunc_1_99.tiff"),
  plot = p,
  width = 6.2,
  height = 4.4,
  dpi = 600,
  compression = "lzw"
)

readr::write_csv(
  weights_plot,
  file.path(TABLES_DIR, "S1R_weight_distribution_trunc_1_99_data.csv")
)

readr::write_csv(
  summary_lines,
  file.path(TABLES_DIR, "S1R_weight_distribution_trunc_1_99_summary.csv")
)

message("Saved:")
message(file.path(FIGURES_DIR, "S1R_weight_distribution_trunc_1_99.png"))
message("==== SUPPLEMENTARY FIGURE S1 COMPLETE ====")