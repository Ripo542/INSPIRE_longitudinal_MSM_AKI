# Script Inventory

## Phase 00: Setup
00_setup_project.R — Initialize paths, load packages, set global options.

## Phase 01: Data Validation
01N_check_inspire_files.R — Verify availability of raw INSPIRE data files.
02N_audit_vitals_resolution.R — Assess temporal resolution and completeness of vital signs.

## Phase 10–14: Cohort Construction
10N_build_base_cohort.R — Define surgical cohort inclusion/exclusion and extract baseline demographics.
11N_build_aki_outcome.R — Retrieve postoperative creatinine and define AKI outcome.
12R_build_baseline_labs.R — Extract preoperative laboratory values.
13R_audit_comorbidities.R — Validate comorbidity coding.
14R_build_comorbidity_flags.R — Code comorbidity indicators and flags.

## Phase 20–24: Longitudinal Data Preparation
20N_build_longitudinal_intervals.R — Create 5-minute interval framework for intraoperative period.
21N_audit_treatment_items.R — Validate fluid and vasopressor ingredient mapping and dosages.
22N_build_treatments_5min.R — Aggregate fluid and vasopressor exposures by interval.
23N_merge_longitudinal_msm_dataset.R — Combine intervals, confounders, treatments, and outcome indicators.
24N_freeze_msm_analysis_dataset.R — Finalize and document the analysis dataset.

## Phase 30–35: MSM Modeling
30N_prepare_msm_modeling_dataset.R — Prepare variables and define lagged confounder structures.
31N_fit_treatment_models_weights.R — Fit treatment allocation models and construct inverse probability weights.
32N_fit_weighted_outcome_model.R — Fit initial weighted logistic regression model.
33N_diagnose_outcome_model_stability.R — Assess model convergence and diagnostics.
34N_fit_final_weighted_outcome_model.R — Fit final weighted outcome model.
34R_fit_revised_weighted_outcome_models.R — Alternative model specifications for sensitivity analysis.
35R_weight_diagnostics.R — Assess weight distributions, truncation effects, and covariate balance.

## Phase 40–53: Results Generation
40*_build_table2_*.R — Compute absolute risk differences.
41*_table1_*.R — Generate baseline characteristics table.
42R_describe_revised_cohorts.R — Cohort summary by subgroup.
43R_surgical_specialty_distribution.R — Surgical specialty distribution.
44R_flowchart_counts.R — Flow diagram documentation.
45R_flowchart_exclusion_breakdown.R — Inclusion/exclusion reasons.
48R_audit_longitudinal_dataset_losses.R — Longitudinal completeness assessment.
49R_table_missing_included_vs_excluded.R — Missingness comparison.
51R_abdominal_restriction_breakdown.R — Subgroup restriction analysis.
52R_vasopressor_type_descriptive.R — Vasopressor use description.
53R_weighted_balance_diagnostics.R — Pre- and post-weighting covariate balance.

## Phase 50+: Sensitivity Analyses
50R_temporal_response_hypotension_treatment.R — Temporal response to hypotension treatment.
S2R_temporal_response_hypotension_treatment.R — Supplementary temporal analyses.

## Figures
32R_forest_primary_model.R — Forest plot for primary model.
32R_forest_primary_weighted_vs_unweighted.R — Comparison of weighted and unweighted estimates.
34N_plot_final_forest.R — Final forest plot.
35N_plot_predicted_risk_exposures.R — Predicted risk curves across exposure gradients.
35R_plot_predicted_risk_primary_abdominal.R — Risk curves for primary abdominal cohort.
F4R_treatment_probability_by_prior_map.R — Treatment probability mapping.
S1N_plot_weight_distribution.R — Weight distribution histogram.
S1R_plot_weight_diagnostics.R — Weight diagnostic plots.
S3R_longitudinal_dag.R — Directed acyclic graph.
