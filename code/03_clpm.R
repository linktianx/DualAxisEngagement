# ======================================================================
# Purpose:
# This script performs cross-lagged panel modeling (CLPM) using lavaan
# to follow up on significant activity–CBCL associations (Pfdr<0.05) identified via LMM.
#
# For each activity–CBCL pair found to be significant in LMM results,
# this script loads corresponding input data and fits a CLPM model to
# examine temporal directional effects between activity and CBCL scores
# across two timepoints, adjusting for baseline covariates.
#
# Output includes:
# - Standardized path estimates with Bonferroni/FDR-adjusted p-values
# - Model fit indices
# - Formatted HTML tables for reporting
# ======================================================================


# -----------------------------
# Setup
# -----------------------------
library(lavaan)
library(kableExtra)
library(dplyr)

# -----------------------------
# Locate all input CLPM data
# -----------------------------
clpm_datapath_list <- list.files(
  path = "/Users/tianx/Documents/academic/Projects/ABCD_lifestyle/activity_cbcl_recheck/CLPM_physical28/results/clpm_results",
  pattern = "input_data.csv",
  recursive = TRUE,
  full.names = TRUE
)

# -----------------------------
# Start loop
# -----------------------------
for (file_path in clpm_datapath_list) {
  cat("\n Processing:", file_path, "\n")

  tryCatch({

    # -----------------------------
    # Parse variable names from path
    # -----------------------------
    path_parts <- strsplit(file_path, "/")[[1]]
    x_var <- path_parts[length(path_parts) - 2]
    last_folder <- path_parts[length(path_parts) - 1]
    last_parts <- strsplit(last_folder, "-")[[1]]
    y_var <- last_parts[1]
    timepoint1 <- last_parts[2]
    timepoint2 <- last_parts[3]

    # -----------------------------
    # Read data and construct model
    # -----------------------------
    df_clpm_input <- read.csv(file_path)
    cv_list <- c("sex", "age", "family_income", "education", "zBMI", "puberty", "race", "site")
    cv_list_time1 <- paste0(cv_list, "_", timepoint1)
    cv_list_time2 <- paste0(cv_list, "_", timepoint2)

    model_desc <- paste0("
      # Cross-lagged paths
      ", x_var, "_", timepoint2, " ~ ", x_var, "_", timepoint1, " + ", y_var, "_", timepoint1, "
      ", y_var, "_", timepoint2, " ~ ", y_var, "_", timepoint1, " + ", x_var, "_", timepoint1, "
      
      # Covariances
      ", x_var, "_", timepoint1, " ~~ ", y_var, "_", timepoint1, "
      ", x_var, "_", timepoint2, " ~~ ", y_var, "_", timepoint2, "
      
      # Covariates for t1 variables
      ", x_var, "_", timepoint1, " ~ ", paste(cv_list_time1, collapse = " + "), "
      ", y_var, "_", timepoint1, " ~ ", paste(cv_list_time1, collapse = " + "), "
      
      # Covariates for t2 variables
      ", x_var, "_", timepoint2, " ~ ", paste(cv_list_time2, collapse = " + "), "
      ", y_var, "_", timepoint2, " ~ ", paste(cv_list_time2, collapse = " + "), "
    ")

    # -----------------------------
    # Run lavaan model
    # -----------------------------
    fit <- sem(model_desc, data = df_clpm_input, missing = "ML")
    fit_summary <- summary(fit, fit.measures = TRUE, standardized = TRUE)

    # -----------------------------
    # Save parameter results
    # -----------------------------
    path_results <- standardizedSolution(fit)
    path_results$p_adjusted_bonferroni <- p.adjust(path_results$pvalue, method = "bonferroni")
    path_results$p_adjusted_fdr <- p.adjust(path_results$pvalue, method = "fdr")

    write.csv(path_results, file.path(dirname(file_path), "lavaan_clpm_results.csv"), row.names = FALSE)

    fit_stats <- as.data.frame(fitMeasures(fit), stringsAsFactors = FALSE)
    write.csv(fit_stats, file.path(dirname(file_path), "lavaan_clpm_fit_eva.csv"), row.names = TRUE)

    # -----------------------------
    # Save HTML table for paths
    # -----------------------------
    path_results_display <- path_results %>%
      filter(op == "~") %>%
      select(lhs, rhs, est.std, pvalue, p_adjusted_bonferroni, p_adjusted_fdr)

    html_table <- path_results_display %>%
      kable(format = "html",
            col.names = c("Dependent", "Predictor", "Standardized β", "p-value", 
                          "Bonferroni-adjusted", "FDR-adjusted")) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)

    save_kable(html_table, file = file.path(dirname(file_path), "path_results_display.html"))

  }, error = function(e) {
    message("Error in file: ", file_path)
    message("   → ", e$message)
  })
}
