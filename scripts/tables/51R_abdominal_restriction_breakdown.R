# =============================================================================
# INSPIRE_HEMO_MSM_AKI | 51R_abdominal_restriction_breakdown.R
# Purpose:
#   Describe the restriction from the longitudinal analytic dataset
#   (n = 44,514 operations) to the elective major abdominal candidate cohort
#   (n = 16,522 operations), without assuming variable names beyond those present
#   in the final 24R longitudinal dataset.
# =============================================================================

source("R/00_setup/00N_setup_project.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

message("==== ABDOMINAL RESTRICTION BREAKDOWN ====")

# -----------------------------------------------------------------------------
# Load final longitudinal analytic dataset
# -----------------------------------------------------------------------------

long <- readRDS(
  file.path(DERIVED_DIR, "24R_revised_longitudinal_msm_with_labs_comorbids.rds")
) %>%
  dplyr::mutate(
    op_id = as.character(op_id),
    subject_id = as.character(subject_id)
  )

message("Rows loaded: ", nrow(long))
message("Operations loaded: ", dplyr::n_distinct(long$op_id))

message("Available columns:")
print(names(long))

# -----------------------------------------------------------------------------
# Operation-level dataset
# -----------------------------------------------------------------------------

op <- long %>%
  dplyr::arrange(op_id, interval_id) %>%
  dplyr::group_by(op_id) %>%
  dplyr::summarise(
    subject_id = dplyr::first(subject_id),
    
    # Final flags already created in 24R
    full_major_noncardiac = dplyr::first(full_major_noncardiac),
    elective_abdominal_major = dplyr::first(elective_abdominal_major),
    revised_primary_abdominal_complete = dplyr::first(revised_primary_abdominal_complete),
    
    # Original/clinical descriptors, if present
    emop = if ("emop" %in% names(long)) dplyr::first(emop) else NA,
    antype = if ("antype" %in% names(long)) dplyr::first(antype) else NA,
    department_clean = if ("department_clean" %in% names(long)) dplyr::first(department_clean) else NA,
    duration_min = if ("duration_min" %in% names(long)) dplyr::first(duration_min) else NA_real_,
    
    baseline_haemoglobin = if ("baseline_haemoglobin" %in% names(long)) dplyr::first(baseline_haemoglobin) else NA_real_,
    baseline_albumin = if ("baseline_albumin" %in% names(long)) dplyr::first(baseline_albumin) else NA_real_,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    included_abdominal_candidate = as.integer(elective_abdominal_major == 1),
    included_primary_abdominal = as.integer(revised_primary_abdominal_complete == 1),
    
    excluded_from_abdominal_candidate = as.integer(elective_abdominal_major != 1),
    
    emop_chr = as.character(emop),
    department_chr = as.character(department_clean),
    antype_chr = as.character(antype),
    
    # Conservative descriptive flags.
    # These are NOT used to redefine the cohort; they are only for describing
    # what is outside elective_abdominal_major.
    non_elective_or_emergency = dplyr::case_when(
      is.na(emop_chr) ~ NA_integer_,
      emop_chr %in% c("1", "TRUE", "true", "T", "Emergency", "emergency") ~ 1L,
      TRUE ~ 0L
    ),
    
    abdominal_specialty_proxy = dplyr::case_when(
      department_chr %in% c("GS", "OG", "UR") ~ 1L,
      is.na(department_chr) ~ NA_integer_,
      TRUE ~ 0L
    )
  )

# -----------------------------------------------------------------------------
# Core counts
# -----------------------------------------------------------------------------

core_counts <- tibble::tibble(
  Step = c(
    "Longitudinal analytic dataset",
    "Elective major abdominal candidate cohort",
    "Primary elective major abdominal cohort",
    "Not in elective major abdominal candidate cohort"
  ),
  Procedures = c(
    dplyr::n_distinct(op$op_id),
    sum(op$included_abdominal_candidate == 1, na.rm = TRUE),
    sum(op$included_primary_abdominal == 1, na.rm = TRUE),
    sum(op$excluded_from_abdominal_candidate == 1, na.rm = TRUE)
  ),
  Patients = c(
    dplyr::n_distinct(op$subject_id),
    dplyr::n_distinct(op$subject_id[op$included_abdominal_candidate == 1]),
    dplyr::n_distinct(op$subject_id[op$included_primary_abdominal == 1]),
    dplyr::n_distinct(op$subject_id[op$excluded_from_abdominal_candidate == 1])
  )
)

# -----------------------------------------------------------------------------
# Mutually exclusive descriptive breakdown of the 44,514 -> 16,522 restriction
# Uses final elective_abdominal_major as source of truth.
# The categories below are descriptive and sequential.
# -----------------------------------------------------------------------------

restriction_breakdown <- op %>%
  dplyr::mutate(
    restriction_category = dplyr::case_when(
      included_abdominal_candidate == 1 ~
        "Included: elective major abdominal candidate cohort",
      
      abdominal_specialty_proxy == 0 ~
        "Excluded from primary cohort definition: non-abdominal surgical specialty",
      
      abdominal_specialty_proxy == 1 & non_elective_or_emergency == 1 ~
        "Excluded from primary cohort definition: abdominal but emergency/non-elective",
      
      abdominal_specialty_proxy == 1 & is.na(non_elective_or_emergency) ~
        "Excluded from primary cohort definition: abdominal specialty with missing urgency flag",
      
      abdominal_specialty_proxy == 1 & non_elective_or_emergency == 0 ~
        "Excluded from primary cohort definition: abdominal specialty but not classified as candidate by final flag",
      
      TRUE ~
        "Excluded from primary cohort definition: unresolved descriptor pattern"
    )
  ) %>%
  dplyr::count(restriction_category, name = "Procedures") %>%
  dplyr::mutate(
    Percent = 100 * Procedures / sum(Procedures)
  ) %>%
  dplyr::arrange(dplyr::desc(Procedures))

# -----------------------------------------------------------------------------
# Surgical specialty distribution inside and outside abdominal candidate cohort
# -----------------------------------------------------------------------------

specialty_distribution <- op %>%
  dplyr::mutate(
    Cohort = dplyr::if_else(
      included_abdominal_candidate == 1,
      "Elective major abdominal candidate cohort",
      "Not in elective major abdominal candidate cohort"
    ),
    Surgical_specialty = dplyr::coalesce(department_chr, "Missing")
  ) %>%
  dplyr::count(Cohort, Surgical_specialty, name = "Procedures") %>%
  dplyr::group_by(Cohort) %>%
  dplyr::mutate(
    Percent = 100 * Procedures / sum(Procedures)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Cohort, dplyr::desc(Procedures))

# -----------------------------------------------------------------------------
# Urgency / emergency descriptor distribution inside and outside candidate cohort
# -----------------------------------------------------------------------------

urgency_distribution <- op %>%
  dplyr::mutate(
    Cohort = dplyr::if_else(
      included_abdominal_candidate == 1,
      "Elective major abdominal candidate cohort",
      "Not in elective major abdominal candidate cohort"
    ),
    Emergency_or_urgent = dplyr::case_when(
      non_elective_or_emergency == 1 ~ "Emergency/non-elective",
      non_elective_or_emergency == 0 ~ "Elective",
      TRUE ~ "Missing"
    )
  ) %>%
  dplyr::count(Cohort, Emergency_or_urgent, name = "Procedures") %>%
  dplyr::group_by(Cohort) %>%
  dplyr::mutate(
    Percent = 100 * Procedures / sum(Procedures)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Cohort, dplyr::desc(Procedures))

# -----------------------------------------------------------------------------
# Final exclusions from abdominal candidate to primary abdominal final cohort
# -----------------------------------------------------------------------------

abdominal_final_exclusions <- op %>%
  dplyr::filter(included_abdominal_candidate == 1) %>%
  dplyr::mutate(
    final_exclusion_reason = dplyr::case_when(
      included_primary_abdominal == 1 ~
        "Included in primary elective major abdominal cohort",
      
      is.na(baseline_haemoglobin) ~
        "Missing baseline hemoglobin",
      
      is.na(baseline_albumin) ~
        "Missing baseline albumin",
      
      TRUE ~
        "Other unresolved final exclusion"
    )
  ) %>%
  dplyr::count(final_exclusion_reason, name = "Procedures") %>%
  dplyr::mutate(
    Percent = 100 * Procedures / sum(Procedures)
  ) %>%
  dplyr::arrange(dplyr::desc(Procedures))

# -----------------------------------------------------------------------------
# Print
# -----------------------------------------------------------------------------

message("==== Core counts ====")
print(tibble::as_tibble(core_counts), n = Inf, width = Inf)

message("==== Restriction breakdown: longitudinal analytic dataset -> abdominal candidate ====")
print(tibble::as_tibble(restriction_breakdown), n = Inf, width = Inf)

message("==== Surgical specialty distribution ====")
print(tibble::as_tibble(specialty_distribution), n = Inf, width = Inf)

message("==== Urgency distribution ====")
print(tibble::as_tibble(urgency_distribution), n = Inf, width = Inf)

message("==== Final abdominal cohort exclusions ====")
print(tibble::as_tibble(abdominal_final_exclusions), n = Inf, width = Inf)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

readr::write_csv(
  core_counts,
  file.path(TABLES_DIR, "51R_abdominal_restriction_core_counts.csv")
)

readr::write_csv(
  restriction_breakdown,
  file.path(TABLES_DIR, "51R_abdominal_restriction_breakdown.csv")
)

readr::write_csv(
  specialty_distribution,
  file.path(TABLES_DIR, "51R_abdominal_restriction_specialty_distribution.csv")
)

readr::write_csv(
  urgency_distribution,
  file.path(TABLES_DIR, "51R_abdominal_restriction_urgency_distribution.csv")
)

readr::write_csv(
  abdominal_final_exclusions,
  file.path(TABLES_DIR, "51R_abdominal_final_exclusions.csv")
)

saveRDS(
  list(
    core_counts = core_counts,
    restriction_breakdown = restriction_breakdown,
    specialty_distribution = specialty_distribution,
    urgency_distribution = urgency_distribution,
    abdominal_final_exclusions = abdominal_final_exclusions
  ),
  file.path(DERIVED_DIR, "51R_abdominal_restriction_breakdown.rds")
)

message("Saved:")
message(file.path(TABLES_DIR, "51R_abdominal_restriction_core_counts.csv"))
message(file.path(TABLES_DIR, "51R_abdominal_restriction_breakdown.csv"))
message(file.path(TABLES_DIR, "51R_abdominal_restriction_specialty_distribution.csv"))
message(file.path(TABLES_DIR, "51R_abdominal_restriction_urgency_distribution.csv"))
message(file.path(TABLES_DIR, "51R_abdominal_final_exclusions.csv"))
message("==== ABDOMINAL RESTRICTION BREAKDOWN COMPLETE ====")