Longitudinal intraoperative hypotension and postoperative acute kidney injury

Marginal structural models for intraoperative hypotension, treatment allocation, and postoperative acute kidney injury.

This repository contains R code used for longitudinal restructuring, inverse probability weighting, and weighted outcome modeling in the INSPIRE perioperative dataset. The analysis addresses time-varying confounding affected by prior treatment using longitudinal marginal structural models.

Contents:

* Longitudinal cohort construction with 5-minute interval restructuring
* Treatment allocation models for fluid and vasopressor administration
* Inverse probability weight construction and diagnostics
* Weighted outcome regression and marginal risk estimation

Data Access

The INSPIRE dataset is not included in this repository. Raw patient-level data are subject to institutional governance and require separate data access agreements.

Software

R version 4.0 or later with packages including tidyverse, data.table, survey, geepack, boot, ggplot2, gridExtra, knitr, and kableExtra. Package versions are documented in logs/sessionInfo.txt.

Configuration

Edit scripts/setup/00_setup_project.R to specify local paths for project_root, data_root, and output directories.

Structure

Scripts are organized into sequential analytical phases numbered 00–53. Directory organization is described in docs/STRUCTURE.md.

Reproducibility

Code and documentation are version-controlled at https://github.com/Ripo542/INSPIRE_longitudinal_MSM_AKI. Package versions are documented in logs/sessionInfo.txt.
