# ======================================================================
# Purpose:
# This script performs linear mixed-effects modeling (LMM) to identify
# associations between activity features (digital and physical) and:
#  (1) CBCL behavioral outcomes
#  (2) brain regional gray matter volumes
#
# Output includes:
# - Standardized LMM coefficients with Wald-based confidence intervals
# - FDR-adjusted p-values for both CBCL and gray matter analyses
# - Separate CSV files for behavioral and neuroimaging outcomes
# ======================================================================


# -----------------------------
# Setup
# -----------------------------
# Load necessary libraries
library(dplyr)
library(readr)
library(stringr)
library(scales)
library(lme4)
library(broom)
library(broom.mixed)

# -----------------------------
# Configurable paths (user-defined)
# -----------------------------
base_dir <- "/Users/tianx/Documents/academic/Projects/ABCD_lifestyle"

# Input data
path_digital_baseline <- file.path(base_dir, "recheck_codes/activity_data/digital/screentime_baseline.csv")
path_physical_baseline <- file.path(base_dir, "recheck_codes/activity_data/physical_ps11_pns17/summary_activity_all_baseline.csv")
path_covariates <- file.path(base_dir, "codes/clean_data/covariates/covariates_all.csv")
path_cbcl <- file.path(base_dir, "codes/clean_data/cbcl/abcd_cbcls01_baseline.csv")
path_grayvol <- file.path(base_dir, "codes/clean_data/dk_fs/baseline/GrayVol_baseline.csv")

# Output directory
res_save_dir <- file.path(base_dir, "activity_cbcl_recheck/code/lmm_r/resutls_physical_28_cbcl11_merge_brain/fdr_bh/lmm_cbcl_baseline")

# Create output directory if not exists
if (!dir.exists(res_save_dir)) {
  dir.create(res_save_dir, recursive = TRUE)
}

# -----------------------------
# Define CBCL variable groups
# -----------------------------
internalizing_problems <- c("cbcl_scr_syn_anxdep_r", "cbcl_scr_syn_withdep_r", "cbcl_scr_syn_somatic_r", "cbcl_scr_syn_internal_r")
externalizing_problems <- c("cbcl_scr_syn_rulebreak_r", "cbcl_scr_syn_aggressive_r", "cbcl_scr_syn_external_r")
problem_behaviors <- c("cbcl_scr_syn_social_r", "cbcl_scr_syn_thought_r", "cbcl_scr_syn_attention_r")
mental_health_problem_total_scale <- c("cbcl_scr_syn_totprob_r")

cbcl_include_list <- c(
  internalizing_problems, 
  externalizing_problems, 
  problem_behaviors, 
  mental_health_problem_total_scale
)

# -----------------------------
# Load and preprocess data
# -----------------------------
df_digital_baseline <- read_csv(path_digital_baseline) %>%
  dplyr::select(src_subject_id, digital_social, digital_nsocial)

df_physical_baseline <- read_csv(path_physical_baseline) %>%
  dplyr::select(src_subject_id, physical_social, physical_nsocial)

df_activity <- inner_join(df_digital_baseline, df_physical_baseline, by = "src_subject_id") %>%
  mutate(src_subject_id = str_replace_all(src_subject_id, "_", ""))

df_covariate <- read_csv(path_covariates) %>%
  mutate(
    src_subject_id = str_replace_all(src_subject_id, "_", ""),
    sex = ifelse(sex == "M", 1, -1)
  )

df_cbcl <- read_csv(path_cbcl) %>%
  dplyr::select(src_subject_id, all_of(cbcl_include_list)) %>%
  mutate(src_subject_id = str_replace_all(src_subject_id, "_", ""))

df_cbcl_std <- df_cbcl
df_cbcl_std[cbcl_include_list] <- scale(df_cbcl[cbcl_include_list])

df_grayvol <- read_csv(path_grayvol)
colnames(df_grayvol)[colnames(df_grayvol) == 'subjectkey'] <- 'src_subject_id'
df_grayvol <- df_grayvol %>%
  mutate(src_subject_id = str_replace_all(src_subject_id, "_", ""))
df_grayvol_std <- df_grayvol
df_grayvol_std[, -1] <- scale(df_grayvol[, -1])

# Standardize covariates
covariate_list <- colnames(df_covariate)[-1]
df_covariate_std <- df_covariate
df_covariate_std[covariate_list] <- scale(df_covariate[covariate_list])


# -----------------------------
# Run linear mixed models : CBCL ~ activity + covariates + random_effect
# -----------------------------
# Model specifications
covariates <- c("sex", "age", "family_income", "education", "zBMI", "puberty", "race")
random_effect <- "(1 | site)"
results_list <- list()
# Loop through each independent variable (activity)
for (independent_var in colnames(df_activity)[-1]) {
  
  df_activity_take <- df_activity %>%
    dplyr::select(src_subject_id, !!sym(independent_var)) %>%
    inner_join(df_cbcl_std, by = "src_subject_id") %>%
    inner_join(df_covariate_std, by = "src_subject_id")

  print(paste("Data shape after merging:", nrow(df_activity_take), "rows"))
  df_activity_take[[independent_var]] <- log10(df_activity_take[[independent_var]] + 1)
  df_activity_take[[independent_var]] <- df_activity_take[[independent_var]] - mean(df_activity_take[[independent_var]], na.rm = TRUE)
  
  for (dv in cbcl_include_list) {
    formula <- as.formula(
      paste(dv, "~", independent_var, "+", paste(covariates, collapse = " + "), "+", random_effect)
    )
    model <- lmer(formula, data = df_activity_take, REML = FALSE)
    coef_table <- broom.mixed::tidy(model, effects = "fixed")
    ind_var_row <- coef_table %>% filter(term == independent_var)
    ci_wald <- confint(model, method = "Wald")
    ci_row <- ci_wald[rownames(ci_wald) == independent_var, , drop = FALSE]
    ind_var_row <- ind_var_row %>%
      mutate(
        independent_var = independent_var,
        dependent_var = dv,
        ci_lower = ci_row[1],
        ci_upper = ci_row[2]
      )
    results_list[[length(results_list) + 1]] <- ind_var_row
  }
}


# Post-processing and save
results_df <- bind_rows(results_list) %>%
  mutate(
    p.value = 2 * pt(-abs(statistic), df = nrow(df_activity_take) - length(covariates) - 1),
    p_value_adjusted = p.adjust(p.value, method = "fdr"),
    significant = p_value_adjusted < 0.05
  )
write_csv(results_df, file.path(res_save_dir, "combined_lmm_CBCL_results_with_ci.csv"))


# -----------------------------
# Run linear mixed models : Greyvol ~ activity + covariates + random_effect
# -----------------------------

covariates <- c("sex", "age", "family_income", "education", "zBMI", "puberty", "race")
random_effect <- "(1 | site)"  
results_list <- list()
for (independent_var in colnames(df_activity)[-1]) {
  # Prepare data
  df_activity_take <- df_activity %>%
    dplyr::select(src_subject_id, !!sym(independent_var)) %>%
    inner_join(df_grayvol_std, by = "src_subject_id") %>%
    inner_join(df_covariate_std, by = "src_subject_id")
  print(paste("Data shape after merging:", nrow(df_activity_take), "rows"))
  df_activity_take[[independent_var]] <- log10(df_activity_take[[independent_var]] + 1)
  df_activity_take[[independent_var]] <- df_activity_take[[independent_var]] - mean(df_activity_take[[independent_var]], na.rm = TRUE)
  dependent_vars <- brain_list
  # Run LMM for each dependent variable
  for (dv in dependent_vars) {
    formula <- as.formula(
      paste(dv, "~", independent_var, "+", paste(covariates, collapse = " + "), "+", random_effect)
    )
    
    model <- lmer(formula, data = df_activity_take, REML = FALSE)
    coef_table <- broom.mixed::tidy(model, effects = "fixed")
    ind_var_row <- coef_table %>% filter(term == independent_var)
    ci_wald <- confint(model, method = "Wald")
    ci_row <- ci_wald[rownames(ci_wald) == independent_var, , drop = FALSE]
    ind_var_row <- ind_var_row %>%
      mutate(
        independent_var = independent_var,
        dependent_var = dv,
        ci_lower = ci_row[1],
        ci_upper = ci_row[2]
      )
    # Append to results
    results_list[[length(results_list) + 1]] <- ind_var_row
  }
}
results_df <- bind_rows(results_list)
results_df <- results_df %>%
  mutate(
    p.value = 2 * pt(-abs(statistic), df = nrow(df_activity_take) - length(covariates) - 1)
  )

results_df <- results_df %>%
  mutate(
    p_value_adjusted = p.adjust(p.value, method = "fdr"),
    significant = p_value_adjusted < 0.05
  )
write.csv(results_df, file = file.path(res_save_dir, "combined_lmm_gw_results_with_ci.csv"), row.names = FALSE)