# ======================================================================

# Purpose:
# This integrated pipeline performs two key stages of analysis:
#   1. Preprocesses raw KSADS diagnostic data at baseline timepoint:
#      - Extracts baseline wave
#      - Replaces special codes
#      - Drops high-missing variables
#      - Retains diagnostic-relevant features
#
#   2. Performs group-level t-tests comparing 4 activity across diagnostic
#      disorder groups defined by KSADS DSM item sets.
#
# Outputs:
#   - Cleaned KSADS baseline dataset
#   - Diagnosis-only filtered dataset
#   - Group-level t-test results with FDR correction
# ======================================================================

# ==================== PATH CONFIGURATION ====================
import os
import re
import pandas as pd
import pingouin as pg
from statsmodels.stats.multitest import multipletests

# Paths (edit as needed)
RAW_DATA_PATH = r"<your_path>\abcd_ksad01.txt"
DICTIONARY_PATH = r"<your_path>\abcd_ksad01.csv"
INTERIM_DATA_DIR = r"<your_path>"
CLEANED_DATA_PATH = r"<your_path>\baseline_cleaned_data.csv"
DROPPED_FEATURES_PATH = r"<your_path>\baseline_dropped_features.csv"
FILTERED_DATA_PATH = r"<your_path>\baseline_diagnosis_data.csv"

# Input for group-level analysis (assumes filtered output reused)
INPUT_CLINICAL_DATA = FILTERED_DATA_PATH
OUTPUT_GROUP_RESULTS = r"<your_path>\group_ttest_results.csv"
OUTPUT_GROUP_EXCLUSIONS = r"<your_path>\group_ttest_exclusions.csv"

# ==================== PREPROCESSING PARAMETERS ====================
TARGET_EVENT = "baseline_year_1_arm_1"

# ==================== DIAGNOSIS GROUP DEFINITIONS ====================
GROUP_DISORDER_DEFINITIONS = {
    # [Same as before, omitted here for brevity]
    'Bipolar': ['ksads_2_837_p', 'ksads_2_835_p', 'ksads_2_836_p'],
    'Suicidality': ['ksads_23_946_p', 'ksads_23_957_p']  # etc...
}

ACTIVITY_VARS = ["digital_social", "digital_nsocial", "physical_social", "physical_nsocial"]

# ==================== STAGE 1: BASELINE PROCESSING ====================
def extract_baseline_data(input_path, output_dir):
    df = pd.read_csv(input_path, delimiter="\t", dtype=str)
    baseline_df = df[df['eventname'] == TARGET_EVENT].copy()
    base_name = os.path.splitext(os.path.basename(input_path))[0]
    output_path = os.path.join(output_dir, f"{base_name}_baseline.txt")
    os.makedirs(output_dir, exist_ok=True)
    baseline_df.to_csv(output_path, sep="\t", index=False)
    return output_path

def clean_special_codes(input_path):
    df = pd.read_csv(input_path, delimiter="\t", dtype=str)
    ksads_cols = [c for c in df.columns if c.startswith('ksads_')]
    df[ksads_cols] = df[ksads_cols].replace({'555', '888', 555, 888, '555.0', '888.0'}, pd.NA)
    return df

def filter_low_variance(df, threshold=0.99):
    drop_cols = df.columns[df.isnull().mean() > threshold]
    return df.drop(columns=drop_cols), list(drop_cols)

def filter_columns_with_dictionary(df, dict_path):
    dict_df = pd.read_csv(dict_path, usecols=['ElementName', 'ElementDescription'])
    diagnosis_cols = dict_df.loc[
        dict_df['ElementDescription'].str.contains('diagnosis', case=False, na=False),
        'ElementName'
    ].tolist()
    core_cols = ['src_subject_id', 'interview_age', 'sex', 'eventname']
    keep = core_cols + [c for c in diagnosis_cols if c in df.columns]
    return df[keep].copy()

def baseline_pipeline():
    raw_path = extract_baseline_data(RAW_DATA_PATH, INTERIM_DATA_DIR)
    df_cleaned = clean_special_codes(raw_path)
    df_final, dropped = filter_low_variance(df_cleaned)
    df_filtered = filter_columns_with_dictionary(df_final, DICTIONARY_PATH)
    df_final.to_csv(CLEANED_DATA_PATH, index=False)
    pd.DataFrame(dropped, columns=['feature']).to_csv(DROPPED_FEATURES_PATH, index=False)
    df_filtered.to_csv(FILTERED_DATA_PATH, index=False)
    return df_filtered

# ==================== STAGE 2: GROUP T-TEST ANALYSIS ====================
def run_group_level_analysis(df, control_mask):
    results, exclusions = [], []
    for disorder, items in GROUP_DISORDER_DEFINITIONS.items():
        case_mask = df[items].any(axis=1)
        for activity in ACTIVITY_VARS:
            ctrl, case = df.loc[control_mask, activity], df.loc[case_mask, activity]
            if len(ctrl) < 2 or len(case) < 2 or ctrl.empty or case.empty:
                exclusions.append({
                    "disorder_group": disorder, "activity_variable": activity,
                    "control_n": len(ctrl), "case_n": len(case),
                    "exclusion_reason": "Too small or empty"
                })
                continue
            stats = pg.ttest(ctrl, case, correction='auto')
            results.append({
                "disorder_group": disorder,
                "activity_variable": activity,
                "t_statistic": stats['T'].values[0],
                "p_value": stats['p-val'].values[0],
                "mean_control": ctrl.mean(),
                "mean_case": case.mean(),
                "control_n": len(ctrl),
                "case_n": len(case),
                "cohens_d": stats['cohen-d'].values[0],
                "statistical_power": stats['power'].values[0]
            })
    return pd.DataFrame(results), exclusions

def apply_multiple_corrections(df):
    out = []
    for var in df['activity_variable'].unique():
        subset = df[df['activity_variable'] == var].copy()
        _, fdr, _, _ = multipletests(subset['p_value'], method='fdr_bh')
        subset['fdr'] = fdr
        out.append(subset)
    return pd.concat(out).reset_index(drop=True)

def save_outputs(results_df, exclusions):
    results_df.to_csv(OUTPUT_GROUP_RESULTS, index=False)
    if exclusions:
        pd.DataFrame(exclusions).to_csv(OUTPUT_GROUP_EXCLUSIONS, index=False)

def group_analysis_pipeline(df_filtered):
    all_items = [i for sublist in GROUP_DISORDER_DEFINITIONS.values() for i in sublist]
    control_mask = (df_filtered[all_items].sum(axis=1) == 0)
    results, exclusions = run_group_level_analysis(df_filtered, control_mask)
    corrected = apply_multiple_corrections(results)
    save_outputs(corrected, exclusions)

# ==================== MASTER EXECUTION ====================
if __name__ == "__main__":
    print("Starting baseline processing...")
    diagnosis_data = baseline_pipeline()
    
    print("Starting group-level t-test analysis...")
    group_analysis_pipeline(diagnosis_data)

    print("All steps complete.")
