# =============================================================================
# 13R_audit_comorbidities.R
# Purpose:
#   Audit whether INSPIRE contains usable comorbidity variables
# =============================================================================

source("R/00_setup/00N_setup_project.R")

library(data.table)
library(dplyr)
library(janitor)

message("==== AUDIT COMORBIDITIES ====")

# -----------------------------------------------------------------------------
# CHECK AVAILABLE FILES
# -----------------------------------------------------------------------------

files <- list.files(RAW_DIR)
print(files)

# -----------------------------------------------------------------------------
# TRY DIAGNOSIS FILE (if exists)
# -----------------------------------------------------------------------------

diag_path <- file.path(RAW_DIR, "diagnosis.csv")

if (!file.exists(diag_path)) {
  message("❌ diagnosis.csv NOT found → comorbidities likely NOT available")
} else {
  
  message("✔ diagnosis.csv found → inspecting")
  
  diag <- fread(diag_path, showProgress = TRUE) %>%
    clean_names()
  
  print(names(diag))
  print(head(diag))
  
  # Identify code column
  possible_code_cols <- names(diag)[grepl("code|icd", names(diag))]
  print(possible_code_cols)
  
  # Basic frequency of codes
  if (length(possible_code_cols) > 0) {
    code_col <- possible_code_cols[1]
    
    code_summary <- diag %>%
      count(.data[[code_col]], sort = TRUE) %>%
      head(50)
    
    message("Top diagnosis codes:")
    print(code_summary)
  }
}

# -----------------------------------------------------------------------------
# CHECK IF COMORBIDITIES ALREADY IN OPERATIONS
# -----------------------------------------------------------------------------

ops_path <- file.path(RAW_DIR, "operations.csv")

ops <- fread(ops_path, showProgress = TRUE) %>%
  clean_names()

message("Operations columns:")
print(names(ops))

# look for obvious comorbidity variables
comorb_cols <- names(ops)[grepl(
  "diab|htn|hypertension|ckd|renal|copd|cardio|stroke",
  names(ops),
  ignore.case = TRUE
)]

message("Detected possible comorbidity columns:")
print(comorb_cols)

# -----------------------------------------------------------------------------
# SUMMARY
# -----------------------------------------------------------------------------

message("==== SUMMARY ====")
message("1) diagnosis.csv present? ", file.exists(diag_path))
message("2) direct comorbidity variables in operations? ", length(comorb_cols) > 0)