# =============================================================================
# 32R_forest_primary_model.R
# =============================================================================

source("R/00_setup/00N_setup_project.R")

library(dplyr)
library(ggplot2)
library(readr)

fit <- readRDS(file.path(DERIVED_DIR, "34R_model_primary_abdominal.rds"))

res <- broom::tidy(fit) %>%
  dplyr::filter(term %in% c(
    "hypotension_time_10min",
    "map_deficit65_per_10",
    "fluid_intervals_per_5",
    "vasopressor_intervals_per_5"
  )) %>%
  dplyr::mutate(
    OR = exp(estimate),
    low = exp(estimate - 1.96 * std.error),
    high = exp(estimate + 1.96 * std.error),
    term = dplyr::recode(
      term,
      hypotension_time_10min = "Hypotension duration (per 10 min)",
      map_deficit65_per_10 = "MAP deficit <65 mmHg (per 10 units)",
      fluid_intervals_per_5 = "Fluid exposure (5-min intervals)",
      vasopressor_intervals_per_5 = "Vasopressor exposure (5-min intervals)"
    )
  )

p <- ggplot(res, aes(x = OR, y = reorder(term, OR))) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey60") +
  geom_errorbarh(aes(xmin = low, xmax = high), height = 0.2, colour = "#1E40AF") +
  geom_point(size = 2.5, colour = "#1E40AF") +
  scale_x_log10() +
  labs(
    x = "Odds ratio (log scale)",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(size = 8.5),
    axis.title.x = element_text(size = 9.5),
    panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3)
  )

print(p)

ggsave(
  file.path(FIGURES_DIR, "32R_forest_primary_model.png"),
  p,
  width = 5,
  height = 3.2,
  dpi = 600
)