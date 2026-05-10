# =============================================================================
# INSPIRE_HEMO_MSM_AKI | S3R_longitudinal_dag.R
# Purpose:
#   Formal dagitty DAG + publication-ready rendering using DAGitty-style roles.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dagitty)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(grid)
})

message("==== LONGITUDINAL DAG ====")

# -----------------------------------------------------------------------------
# Formal DAG object
# -----------------------------------------------------------------------------

dag_code <- "
dag {
B
H_prev
T_prev
T_current
H_current
AKI
U

B -> H_prev
B -> T_prev
B -> T_current
B -> AKI

T_prev -> H_prev
H_prev -> T_current
T_current -> H_current

H_prev -> AKI
H_current -> AKI
T_current -> AKI

U -> H_prev
U -> T_current
U -> AKI
}
"

dag <- dagitty::dagitty(dag_code)

message("DAG object created successfully.")
print(dag)

# -----------------------------------------------------------------------------
# Nodes
# -----------------------------------------------------------------------------

nodes <- tibble::tribble(
  ~name,        ~x, ~y,   ~label,                                                                 ~role,
  "B",          0,  2.4, "Baseline covariates\nage, sex, ASA, BMI,\ncreatinine, hemoglobin, albumin,\ncomorbidities, specialty", "ancestor_both",
  "T_prev",     2,  1.1, "Prior treatment\nfluids / vasopressors\nbefore interval t",              "ancestor_both",
  "H_prev",     2,  3.7, "Prior MAP state\nMAP / hypotension\nat t−1",                              "ancestor_both",
  "T_current",  4,  2.4, "Current treatment\nfluids / vasopressors\nat interval t",                 "exposure",
  "H_current",  6,  2.4, "Updated MAP state\nMAP / hypotension\nat interval t",                     "other",
  "AKI",        8,  2.4, "Postoperative AKI",                                                       "outcome",
  "U",          4,  0.1, "Unmeasured clinical context\n/ treatment indication",                     "unobserved"
)
# -----------------------------------------------------------------------------
# Edges
# -----------------------------------------------------------------------------

edges <- tibble::tribble(
  ~from,        ~to,          ~curvature,
  "B",         "H_prev",      0.08,
  "B",         "T_prev",     -0.08,
  "B",         "T_current",   0.00,
  "B",         "AKI",         0.30,
  
  "T_prev",    "H_prev",      0.15,
  "H_prev",    "T_current",   0.00,
  "T_current", "H_current",   0.00,
  
  "H_prev",    "AKI",         0.18,
  "H_current", "AKI",         0.00,
  "T_current", "AKI",        -0.16,
  
  "U",         "H_prev",      0.18,
  "U",         "T_current",   0.00,
  "U",         "AKI",        -0.24
) %>%
  dplyr::left_join(
    nodes %>% dplyr::select(from = name, x_from = x, y_from = y),
    by = "from"
  ) %>%
  dplyr::left_join(
    nodes %>% dplyr::select(to = name, x_to = x, y_to = y),
    by = "to"
  )

# -----------------------------------------------------------------------------
# DAGitty-style colours, softened for publication
# -----------------------------------------------------------------------------

fill_values <- c(
  exposure      = "#EEF46B",
  outcome       = "#7ED4E6",
  ancestor_both = "#F7C7CC",
  unobserved    = "#F1F5F9",
  other         = "#DDE3EA"
)

edges_plot <- edges %>%
  dplyr::filter(
    is.finite(x_from),
    is.finite(y_from),
    is.finite(x_to),
    is.finite(y_to),
    !(x_from == x_to & y_from == y_to)
  ) %>%
  dplyr::mutate(
    edge_type = dplyr::if_else(from == "U", "unmeasured", "measured")
  )

nodes_regular <- nodes %>% dplyr::filter(name != "AKI")
nodes_aki <- nodes %>% dplyr::filter(name == "AKI")

p <- ggplot() +
  geom_curve(
    data = edges_plot %>% dplyr::filter(edge_type == "measured"),
    aes(
      x = x_from,
      y = y_from,
      xend = x_to,
      yend = y_to,
      curvature = curvature
    ),
    colour = "grey38",
    linewidth = 0.44,
    arrow = grid::arrow(
      length = grid::unit(5.2, "pt"),
      type = "closed"
    ),
    lineend = "round"
  ) +
  geom_curve(
    data = edges_plot %>% dplyr::filter(edge_type == "unmeasured"),
    aes(
      x = x_from,
      y = y_from,
      xend = x_to,
      yend = y_to,
      curvature = curvature
    ),
    colour = "grey48",
    linewidth = 0.40,
    linetype = "dashed",
    arrow = grid::arrow(
      length = grid::unit(5.2, "pt"),
      type = "closed"
    ),
    lineend = "round"
  ) +
  geom_label(
    data = nodes_regular,
    aes(
      x = x,
      y = y,
      label = label,
      fill = role
    ),
    colour = "black",
    size = 2.55,
    lineheight = 0.88,
    label.size = 0,
    label.padding = grid::unit(0.34, "lines"),
    label.r = grid::unit(0.18, "lines")
  ) +
  geom_label(
    data = nodes_aki,
    aes(
      x = x,
      y = y,
      label = label,
      fill = role
    ),
    colour = "black",
    size = 2.8,
    lineheight = 0.88,
    label.size = 0,
    label.padding = grid::unit(0.42, "lines"),
    label.r = grid::unit(0.20, "lines")
  ) +
  scale_fill_manual(values = fill_values) +
  coord_cartesian(
    xlim = c(-1.05, 9.05),
    ylim = c(-0.70, 4.45),
    clip = "off"
  ) +
  theme_void(base_size = 10) +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 14, 10, 14)
  )

print(p)
# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

ggsave(
  filename = file.path(FIGURES_DIR, "S3R_longitudinal_dag.png"),
  plot = p,
  width = 9.2,
  height = 5.0,
  dpi = 600
)

ggsave(
  filename = file.path(FIGURES_DIR, "S3R_longitudinal_dag.pdf"),
  plot = p,
  width = 9.2,
  height = 5.0
)

ggsave(
  filename = file.path(FIGURES_DIR, "S3R_longitudinal_dag.tiff"),
  plot = p,
  width = 9.2,
  height = 5.0,
  dpi = 600,
  compression = "lzw"
)

writeLines(
  dag_code,
  con = file.path(QC_DIR, "S3R_longitudinal_dag_dagitty_code.txt")
)

readr::write_csv(
  nodes,
  file.path(QC_DIR, "S3R_longitudinal_dag_nodes.csv")
)

readr::write_csv(
  edges,
  file.path(QC_DIR, "S3R_longitudinal_dag_edges.csv")
)

message("Saved:")
message(file.path(FIGURES_DIR, "S3R_longitudinal_dag.png"))
message(file.path(QC_DIR, "S3R_longitudinal_dag_dagitty_code.txt"))
message("==== LONGITUDINAL DAG COMPLETE ====")