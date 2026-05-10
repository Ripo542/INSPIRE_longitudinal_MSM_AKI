# =============================================================================
# Plot weight distributions (operation-level)
# =============================================================================

source("R/00_setup/00N_setup_project.R")

opw <- readRDS(file.path(DERIVED_DIR, "31N_operation_final_weights.rds"))

# -----------------------------------------------------------------------------
# Prepare data (log-scale for untruncated)
# -----------------------------------------------------------------------------

library(ggplot2)

# avoid log issues
opw <- opw %>%
  dplyr::mutate(
    final_sw_joint_log = log10(final_sw_joint + 1e-10)
  )

# -----------------------------------------------------------------------------
# PLOT 1 — Untruncated (log scale)
# -----------------------------------------------------------------------------

p1 <- ggplot(opw, aes(x = final_sw_joint_log)) +
  geom_histogram(bins = 100) +
  labs(
    title = "Distribution of untruncated stabilized weights (log scale)",
    x = "log10(weight)",
    y = "Count"
  ) +
  theme_minimal()

# -----------------------------------------------------------------------------
# PLOT 2 — Truncated 1–99
# -----------------------------------------------------------------------------

p2 <- ggplot(opw, aes(x = final_sw_joint_trunc_1_99)) +
  geom_histogram(bins = 100) +
  labs(
    title = "Distribution of stabilized weights (1st–99th percentile truncation)",
    x = "Weight",
    y = "Count"
  ) +
  theme_minimal()

# -----------------------------------------------------------------------------
# SAVE
# -----------------------------------------------------------------------------

ggsave(
  filename = file.path(FIGURES_DIR, "weights_untruncated_log.png"),
  plot = p1,
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(FIGURES_DIR, "weights_truncated_1_99.png"),
  plot = p2,
  width = 6,
  height = 4,
  dpi = 300
)

message("Saved figures:")
message(file.path(FIGURES_DIR, "weights_untruncated_log.png"))
message(file.path(FIGURES_DIR, "weights_truncated_1_99.png"))