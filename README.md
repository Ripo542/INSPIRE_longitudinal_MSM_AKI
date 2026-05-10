# Longitudinal intraoperative hemodynamics and acute kidney injury

Marginal structural models for intraoperative hemodynamic exposures and postoperative acute kidney injury.

This repository contains R analytical code implementing inverse probability weighting to examine associations between intraoperative fluid and vasopressor administration and postoperative acute kidney injury in the INSPIRE perioperative dataset. The analysis addresses time-varying confounding affected by prior treatment using longitudinal marginal structural models.

Contents:
- Longitudinal cohort construction with 5-minute interval restructuring
- Treatment allocation models for fluid and vasopressor administration
- Inverse probability weight construction and diagnostics
- Weighted outcome regression and marginal risk estimation

## Data Access

The INSPIRE dataset is not included in this repository. Raw patient-level data are subject to institutional governance and require separate data access agreement. This repository contains analytical code.

## Software

R version 4.0 or later with packages: tidyverse, data.table, survey, geepack, boot, ggplot2, gridExtra, knitr, kableExtra. Package versions documented in logs/sessionInfo.txt.

## Configuration

Edit scripts/setup/00_setup_project.R to specify local paths for project_root, data_root, and output directories.

## Structure

Scripts are organized into phases numbered 00–53, executed in sequential order. See docs/STRUCTURE.md for directory organization.

## Outputs

Tabular results in outputs/tables/, figures in outputs/figures/.

## Reproducibility

All code is version-controlled at https://github.com/Ripo542/INSPIRE_longitudinal_MSM_AKI Package versions are documented in logs/sessionInfo.txt.
