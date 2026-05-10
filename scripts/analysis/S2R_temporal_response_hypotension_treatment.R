# =============================================================================
# INSPIRE_HEMO_MSM_AKI | S2R_temporal_response_hypotension_treatment.R
# Purpose:
#   Supplementary Figure S2:
#   Treatment probability according to hypotension in the previous interval.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(scales)
})

message("==== SUPPLEMENTARY FIGURE S2: TEMPORAL RESPONSE ====")

summary_tbl <- readr::read_csv(
  file.path(TABLES_DIR, "50R_temporal_response_summary.csv"),
  show_col_types = FALSE
)

plot_dat <- summary_tbl %>%
  dplyr::mutate(
    Metric = dplyr::recode(
      Metric,
      "Fluid administration" = "Fluid\nadministration",
      "Vasopressor administration" = "Vasopressor\nadministration",
      "Combined fluid + vasopressor" = "Combined fluid\n+ vasopressor"
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(Hypotension_prev, No_hypotension_prev),
    names_to = "Previous_state",
    values_to = "Probability"
  ) %>%
  dplyr::mutate(
    Previous_state = dplyr::recode(
      Previous_state,
      "Hypotension_prev" = "Previous hypotension",
      "No_hypotension_prev" = "No previous hypotension"
    ),
    Previous_state = factor(
      Previous_state,
      levels = c("No previous hypotension", "Previous hypotension")
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Fluid\nadministration",
        "Vasopressor\nadministration",
        "Combined fluid\n+ vasopressor"
      )
    )
  )

segment_dat <- plot_dat %>%
  tidyr::pivot_wider(
    names_from = Previous_state,
    values_from = Probability
  )

p <- ggplot(
  plot_dat,
  aes(
    x = Probability,
    y = Previous_state,
    fill = Previous_state
  )
) +
  geom_col(
    width = 0.55,
    colour = NA
  ) +
  facet_wrap(
    ~ Metric,
    ncol = 1,
    scales = "free_x",
    strip.position = "top"
  ) +
  scale_fill_manual(
    values = c(
      "No previous hypotension" = "#D8E2EE",
      "Previous hypotension" = "#1E40AF"
    )
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 0.1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Probability of treatment in the current 5-min interval",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "plain", margin = margin(b = 4)),
    axis.text = element_text(size = 8.5, colour = "black"),
    axis.title.x = element_text(size = 9.5, margin = margin(t = 7)),
    axis.line = element_line(colour = "black", linewidth = 0.32),
    axis.ticks = element_line(colour = "black", linewidth = 0.28),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.28),
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 5)
  )

ggsave(
  filename = file.path(FIGURES_DIR, "S2R_temporal_response_hypotension_treatment.png"),
  plot = p,
  width = 5.8,
  height = 3.4,
  dpi = 600
)

ggsave(
  filename = file.path(FIGURES_DIR, "S2R_temporal_response_hypotension_treatment.pdf"),
  plot = p,
  width = 5.8,
  height = 3.4
)

ggsave(
  filename = file.path(FIGURES_DIR, "S2R_temporal_response_hypotension_treatment.tiff"),
  plot = p,
  width = 5.8,
  height = 3.4,
  dpi = 600,
  compression = "lzw"
)

readr::write_csv(
  plot_dat,
  file.path(TABLES_DIR, "S2R_temporal_response_hypotension_treatment_plot_data.csv")
)

message("Saved:")
message(file.path(FIGURES_DIR, "S2R_temporal_response_hypotension_treatment.png"))
message("==== SUPPLEMENTARY FIGURE S2 COMPLETE ====")