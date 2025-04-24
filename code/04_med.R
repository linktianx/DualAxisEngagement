# ======================================================================
# Purpose:
# This script conducts mediation analysis using digital nonsocial activity
# (digital_nsocial_baseline) as the independent variable. It follows up
# on significant activity–CBCL associations identified through CLPM.
#
# The mediator is the total gray matter volume (GMV) across significant
# brain regions identified through LMM. For each selected CBCL outcome,
# the mediation model estimates:
#   - Path a: IV (activity) → Mediator (brain GMV)
#   - Path b: Mediator → DV (CBCL outcome)
#   - Path c: Total effect of IV on DV
#   - Path c′: Direct effect
#   - Indirect effect (a × b), tested via bootstrap with 5,000 simulations
#
# Outputs include:
# - Text summaries for paths a, b, c, c′, and indirect effect
# - Full lavaan model summaries
# - Confidence intervals and p-values for indirect effect
#
# This script uses digital nonsocial activity as a worked example,
# ======================================================================

# -----------------------------
# Load libraries
# -----------------------------
library(mediation)
library(readr)

# -----------------------------
# Define data and analysis setup
# -----------------------------
input_data_dir <- "/Users/tianx/Documents/academic/Projects/ABCD_lifestyle/activity_cbcl_recheck/mediation_PS28_rlmm/fs_baseline/DNS_GrayVol_CBCL"
df_tot <- read_csv(file.path(input_data_dir, "input_data.csv"))

treatment <- "digital_nsocial_baseline"  # Independent variable
mediator <- "gvsum"                    # Mediator from LMM-based results
time_label <- "_f1yr"

# Outcome CBCL variables
outcome_list_base <- c(
  "cbcl_scr_syn_rulebreak_r",
  "cbcl_scr_syn_attention_r",
  "cbcl_scr_syn_totprob_r",
  "cbcl_scr_syn_thought_r",
  "cbcl_scr_syn_social_r"
)
outcome_list <- paste0(outcome_list_base, time_label)
# Covariates (baseline)
covariates_list <- c("sex", "age", "family_income", "education", "zBMI", "puberty", "race", "site")
covariates <- paste0(covariates_list, "_baseline")

# -----------------------------
# Run mediation analysis per outcome
# -----------------------------
for (outcome in outcome_list) {
  
  # Create save folder
  res_save_dir <- file.path(input_data_dir, "med_res", outcome)
  if (!dir.exists(res_save_dir)) dir.create(res_save_dir, recursive = TRUE)
  
  # Path a: IV → Mediator
  mediator_formula <- as.formula(
    paste(mediator, "~", treatment, "+", paste(covariates, collapse = " + "))
  )
  model_a <- lm(mediator_formula, data = df_tot)
  summary_a <- summary(model_a)
  path_a <- coef(model_a)[treatment]
  path_a_pval <- coef(summary_a)[treatment, "Pr(>|t|)"]
  writeLines(capture.output(summary_a), file.path(res_save_dir, "summary_model_a.txt"))
  
  # Path b & c': Mediator + IV → DV
  outcome_formula <- as.formula(
    paste(outcome, "~", treatment, "+", mediator, "+", paste(covariates, collapse = " + "))
  )
  model_b_cprime <- lm(outcome_formula, data = df_tot)
  summary_b_cprime <- summary(model_b_cprime)
  path_b <- coef(model_b_cprime)[mediator]
  path_b_pval <- coef(summary_b_cprime)[mediator, "Pr(>|t|)"]
  path_c_prime <- coef(model_b_cprime)[treatment]
  path_c_prime_pval <- coef(summary_b_cprime)[treatment, "Pr(>|t|)"]
  writeLines(capture.output(summary_b_cprime), file.path(res_save_dir, "summary_model_b_cprime.txt"))
  
  # Path c: Total effect (IV → DV)
  total_effect_formula <- as.formula(
    paste(outcome, "~", treatment, "+", paste(covariates, collapse = " + "))
  )
  model_c <- lm(total_effect_formula, data = df_tot)
  summary_c <- summary(model_c)
  path_c <- coef(model_c)[treatment]
  path_c_pval <- coef(summary_c)[treatment, "Pr(>|t|)"]
  writeLines(capture.output(summary_c), file.path(res_save_dir, "summary_model_c.txt"))
  
  # Mediation test (a × b) via bootstrap
  medM1_result <- mediate(
    model.m = model_a,
    model.y = model_b_cprime,
    treat = treatment,
    mediator = mediator,
    sims = 5000,
    boot = TRUE,
    boot.ci.type = "bca"
  )
  summary_med <- summary(medM1_result)
  writeLines(capture.output(summary_med), file.path(res_save_dir, "summary_mediation.txt"))
  
  # Extract mediation effects
  indirect_effect <- summary_med$d.avg
  indirect_ci <- summary_med$d.avg.ci
  indirect_pval <- summary_med$d.avg.p
  
  # Create text summary
  effect_summary <- paste(
    "Path a (", treatment, " → ", mediator, "): ", round(path_a, 4), ", p = ", format.pval(path_a_pval), "\n",
    "Path b (", mediator, " → ", outcome, "): ", round(path_b, 4), ", p = ", format.pval(path_b_pval), "\n",
    "Path c' (direct): ", round(path_c_prime, 4), ", p = ", format.pval(path_c_prime_pval), "\n",
    "Path c (total): ", round(path_c, 4), ", p = ", format.pval(path_c_pval), "\n",
    "Indirect effect (a × b): ", round(indirect_effect, 4),
    " [95% CI: ", round(indirect_ci[1], 4), ", ", round(indirect_ci[2], 4), "], p = ", format.pval(indirect_pval), "\n",
    sep = ""
  )
  
  writeLines(effect_summary, file.path(res_save_dir, "path_effects_summary.txt"))
  cat("Completed mediation for outcome:", outcome, "\n")
}
