# utils.py
import pandas as pd
import numpy as np
import json
import warnings

def clean_column_names(columns):
    return columns.str.strip().str.replace('"', '', regex=False).str.replace("'", '', regex=False)

def clean_cell_values(series):
    series = series.astype(object)
    series = series.where(series.notna(), '')
    series = series.astype(str)
    try:
        return series.replace(r'^\s*["\']{0,2}\s*$', np.nan, regex=True).str.strip(' "\'')
    except AttributeError:
        return series

def robust_load_csv(filepath, expected_columns=None):
    try:
        df = pd.read_csv(
            filepath, encoding="ISO-8859-1", engine="python",
            sep=",", quotechar='"', on_bad_lines="skip"
        )
    except Exception as e:
        raise RuntimeError(f"[ERROR] Could not load {filepath}: {e}")

    df.columns = clean_column_names(df.columns)
    df = df.loc[:, ~df.columns.duplicated()]

    for col in df.columns:
        df[col] = clean_cell_values(df[col])

    if expected_columns and not expected_columns.issubset(set(df.columns)):
        missing = expected_columns - set(df.columns)
        raise ValueError(f"[ERROR] Missing expected columns: {missing}")

    return df

def apply_fallback_names(df):
    for col in ['Name', 'Formula', 'Calc. MW']:
        if col not in df.columns:
            df[col] = np.nan
        else:
            df[col] = clean_cell_values(df[col])
    df['Name'] = df['Name'].fillna(df['Formula']).fillna(df['Calc. MW'])
    return df

def ensure_column(df, col, default=''):
    if col not in df.columns:
        df[col] = default
    else:
        df[col] = clean_cell_values(df[col])
    return df

def load_and_prepare_data():
    # Load data and prepare for visualization
    distance_df = robust_load_csv('by_distance_named.csv', expected_columns={'Compounds ID', 'Calc. MW', 'Name'})
    distance_df = apply_fallback_names(distance_df)
    gold_ids = distance_df['Compounds ID'].dropna().astype(str).tolist()

    raw_data = robust_load_csv('ReSpleen.csv')
    raw_data = apply_fallback_names(raw_data)
    raw_data = ensure_column(raw_data, 'm/z')
    raw_data = ensure_column(raw_data, 'RT [min]')
    raw_data['Compounds ID'] = raw_data['Compounds ID'].astype(str)
    raw_data['Gold'] = raw_data['Compounds ID'].isin(gold_ids)

    with open("column_mapping.json") as f:
        mapping_data = json.load(f)
    comparisons = mapping_data.get("ReSpleen.csv", [])
    
    if not comparisons:
        raise ValueError("No comparisons found in column_mapping.json")

    # Pre-compute data for all comparisons
    for entry in comparisons:
        fc_col = entry.get("fold_change_col")
        pv_col = entry.get("p_value_col")

        if fc_col in raw_data.columns and pv_col in raw_data.columns:
            raw_data[fc_col] = pd.to_numeric(clean_cell_values(raw_data[fc_col]), errors='coerce')
            raw_data[pv_col] = pd.to_numeric(clean_cell_values(raw_data[pv_col]), errors='coerce')
            raw_data[f'-Log10({pv_col})'] = -np.log10(raw_data[pv_col])
            raw_data[f'{fc_col}_sig_up'] = (raw_data[fc_col] > 0.5) & (raw_data[pv_col] < 0.05)
            raw_data[f'{fc_col}_sig_down'] = (raw_data[fc_col] < -0.5) & (raw_data[pv_col] < 0.05)

    return raw_data, comparisons
