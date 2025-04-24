# DualAxisEngagement

This repository contains the analysis pipeline for the manuscript:
> **Physical-Digital and Social-Nonsocial Extracurricular Engagement: Differential Effects on Brain Development and Psychological Outcomes in Children**


## Code Overview (located in `code/`)
| Script | Description |
|--------|-------------|
| `01_lmm.R` | Fits linear mixed models to estimate associations between activity engagement and CBCL/brain outcomes at baseline |
| `02_diag_ttest.py` | Performs group-level t-tests comparing diagnostic subgroups (DSM-5 based) on different activity variables |
| `03_clpm.R` | Implements cross-lagged panel modeling (CLPM) for significant activity–CBCL pairs to assess directionality |
| `04_med.R` | Conducts mediation analysis to test whether gray matter volume mediates the effect of activity on CBCL outcomes |


Paths and variable names should be adapted based on user-specific ABCD datasets.
