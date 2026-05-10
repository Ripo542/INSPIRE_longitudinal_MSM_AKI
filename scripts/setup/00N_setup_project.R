
# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 00N_setup_project.R
# =============================================================================

PROJECT_NAME <- "INSPIRE_HEMO_MSM_AKI"

find_project_root <- function(start = getwd(), marker = "README.md") {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, marker)) && basename(path) == PROJECT_NAME) return(path)
    parent <- dirname(path)
    if (identical(parent, path)) stop("Project root not found. Open/set working directory inside: ", PROJECT_NAME)
    path <- parent
  }
}

ROOT_DIR <- find_project_root()

DATA_DIR      <- file.path(ROOT_DIR, "data")
LOCAL_RAW_DIR <- file.path(DATA_DIR, "raw")
DERIVED_DIR   <- file.path(DATA_DIR, "derived")
FROZEN_DIR    <- file.path(DATA_DIR, "frozen")

R_DIR        <- file.path(ROOT_DIR, "R")
OUTPUTS_DIR  <- file.path(ROOT_DIR, "outputs")
QC_DIR       <- file.path(OUTPUTS_DIR, "qc")
TABLES_DIR   <- file.path(OUTPUTS_DIR, "tables")
FIGURES_DIR  <- file.path(OUTPUTS_DIR, "figures")
LOGS_DIR     <- file.path(OUTPUTS_DIR, "logs")

INSPIRE_SOURCE_DIR <- "/Users/javi/Documents/INSPIRE/inspire-a-publicly-available-research-dataset-for-perioperative-medicine-1.3"

if (dir.exists(INSPIRE_SOURCE_DIR)) {
  RAW_DIR <- INSPIRE_SOURCE_DIR
} else {
  RAW_DIR <- LOCAL_RAW_DIR
  warning("External INSPIRE_SOURCE_DIR not found. Falling back to local RAW_DIR: ", RAW_DIR)
}

dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FROZEN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOGS_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(20260505)

options(
  stringsAsFactors = FALSE,
  scipen = 999,
  width = 120
)

required_packages <- c(
  "data.table", "dplyr", "tidyr", "readr", "arrow", "jsonlite",
  "lubridate", "janitor", "ggplot2", "broom", "sandwich",
  "lmtest", "cobalt", "WeightIt", "purrr", "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the pipeline."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(arrow)
  library(jsonlite)
  library(lubridate)
  library(janitor)
  library(ggplot2)
  library(purrr)
  library(tibble)
})

start_log <- function(script_name) {
  log_file <- file.path(
    LOGS_DIR,
    paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", script_name, ".log")
  )
  sink(log_file, split = TRUE)
  sink(log_file, type = "message", append = TRUE)
  message("==== START: ", script_name, " ====")
  message("Time: ", Sys.time())
  message("Root: ", ROOT_DIR)
  message("Raw data: ", RAW_DIR)
  message("R version: ", R.version.string)
  return(log_file)
}

end_log <- function(log_file = NULL) {
  message("==== END ====")
  message("Time: ", Sys.time())
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  invisible(log_file)
}

write_qc_json <- function(x, filename) {
  jsonlite::write_json(
    x,
    path = file.path(QC_DIR, filename),
    pretty = TRUE,
    auto_unbox = TRUE,
    null = "null"
  )
}

message("Setup loaded successfully.")
message("ROOT_DIR: ", ROOT_DIR)
message("RAW_DIR:  ", RAW_DIR)

