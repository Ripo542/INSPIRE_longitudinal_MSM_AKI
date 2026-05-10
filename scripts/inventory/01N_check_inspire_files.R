# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 01N_check_inspire_files.R
# Purpose:
#   Inventory raw INSPIRE files and verify availability + minimal structure.
# =============================================================================

SCRIPT <- "01N_check_inspire_files"

source("R/00_setup/00N_setup_project.R")

message("==== INSPIRE RAW FILE INVENTORY ====")

expected_files <- c(
  "operations.csv",
  "vitals.csv",
  "labs.csv",
  "medications.csv"
)

paths <- file.path(RAW_DIR, expected_files)

file_inventory <- tibble::tibble(
  file   = expected_files,
  path   = paths,
  exists = file.exists(paths),
  size_mb = ifelse(
    file.exists(paths),
    round(file.info(paths)$size / 1024^2, 2),
    NA_real_
  )
)

print(file_inventory)

missing_files <- file_inventory$file[!file_inventory$exists]

if (length(missing_files) > 0) {
  stop(
    "Missing required INSPIRE files:\n",
    paste(missing_files, collapse = "\n")
  )
}

read_header <- function(path) {
  names(readr::read_csv(path, n_max = 0, show_col_types = FALSE))
}

headers <- list(
  operations  = read_header(file.path(RAW_DIR, "operations.csv")),
  vitals      = read_header(file.path(RAW_DIR, "vitals.csv")),
  labs        = read_header(file.path(RAW_DIR, "labs.csv")),
  medications = read_header(file.path(RAW_DIR, "medications.csv"))
)

required_columns <- list(
  operations  = c("op_id", "subject_id"),
  vitals      = c("op_id", "subject_id", "chart_time", "item_name", "value"),
  labs        = c("subject_id", "chart_time", "item_name", "value"),
  medications = c("subject_id", "chart_time", "drug_name")
)
column_audit <- dplyr::bind_rows(lapply(names(required_columns), function(dataset_name) {
  cols <- required_columns[[dataset_name]]
  available_cols <- headers[[dataset_name]]
  
  tibble::tibble(
    dataset = dataset_name,
    required_column = cols,
    present = cols %in% available_cols
  )
}))

print(column_audit)

missing_cols <- column_audit %>% dplyr::filter(!present)

if (nrow(missing_cols) > 0) {
  stop(
    "Missing required columns:\n",
    paste(
      paste0(missing_cols$dataset, "::", missing_cols$required_column),
      collapse = "\n"
    )
  )
}

qc <- list(
  script = SCRIPT,
  timestamp = as.character(Sys.time()),
  raw_dir = RAW_DIR,
  file_inventory = file_inventory,
  headers = headers,
  column_audit = column_audit
)

write_qc_json(qc, "01N_check_inspire_files.json")

message("==== INVENTORY OK ====")