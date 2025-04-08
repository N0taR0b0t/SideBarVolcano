import pandas as pd
import numpy as np
import json
from collections import defaultdict

def clean_column_names(columns):
    """Standardizes and strips unnecessary characters from column names."""
    return columns.str.strip().str.replace('"', '', regex=False).str.replace("'", '', regex=False)

def clean_numeric_column(series, col_name):
    """Attempts to clean and convert a column to numeric values."""
    cleaned = series.astype(str).str.strip().str.strip('"').str.strip()
    numeric = pd.to_numeric(cleaned, errors='coerce')
    num_non_numeric = numeric.isna().sum()
    if num_non_numeric == len(series):
        print(f"\nWARNING: ALL VALUES in column '{col_name}' failed numeric conversion!")
        sample_values = series.dropna().unique()[:5]
        print(f"Sample values: {sample_values}")
    elif num_non_numeric > 0:
        print(f"\nWARNING: {num_non_numeric} non-numeric entries detected in column '{col_name}'")
        print(f"Sample values: {series[numeric.isna()].dropna().unique()[:5]}")
    return numeric

def strip_nested_quotes(val):
    """Strips layers of surrounding quotes and whitespace. Converts empty results to NaN."""
    if pd.isna(val):
        return np.nan
    val = str(val).strip()
    while val.startswith(("'", '"')) and val.endswith(("'", '"')) and len(val) >= 2:
        val = val[1:-1].strip()
    return val if val else np.nan

def compute_compound_distances(file: str, mapping_file: str = "column_mapping.json",
                               max_results: int = 200, distance_threshold: float = 7.5) -> pd.DataFrame:
    with open(mapping_file, "r") as f:
        column_mapping = json.load(f)

    try:
        df = pd.read_csv(file, encoding="ISO-8859-1", engine="python", sep=",", quotechar='"', on_bad_lines="skip")
    except FileNotFoundError:
        raise RuntimeError(f"File '{file}' not found. Please ensure it is in the working directory.")
    except Exception as e:
        raise RuntimeError(f"Error reading file '{file}': {e}")

    df.columns = clean_column_names(df.columns)
    print(f"\n[DEBUG] Number of columns read from '{file}': {len(df.columns)}")
    print("\nColumns read from CSV:")
    for col in df.columns:
        print(f"  - {col}")

    comparisons = []
    for entry in column_mapping.get(file, []):
        fc = entry["fold_change_col"].strip().replace('"', '').replace("'", "")
        pv = entry["p_value_col"].strip().replace('"', '').replace("'", "")
        comparisons.append((fc, pv))

    print("\nComparisons extracted from mapping:")
    for fc, pv in comparisons:
        print(f"  - Fold Change: {fc}, P-value: {pv}")

    compound_distances = defaultdict(list)

    for log2fc_col, pval_col in comparisons:
        print(f"\n[DEBUG] Processing comparison: {log2fc_col} vs {pval_col}")

        if log2fc_col not in df.columns or pval_col not in df.columns:
            print(f"[WARNING] Skipping comparison: Missing column(s): '{log2fc_col}' or '{pval_col}' in '{file}'")
            continue

        df[pval_col] = clean_numeric_column(df[pval_col], pval_col)
        df[log2fc_col] = clean_numeric_column(df[log2fc_col], log2fc_col)
        df['-log10(p-value)'] = -np.log10(df[pval_col])

        valid_df = df.dropna(subset=[log2fc_col, '-log10(p-value)', 'Compounds ID'])

        if valid_df.empty:
            print(f"[WARNING] No valid data found for comparison '{log2fc_col}' vs '{pval_col}'. Skipping.")
            continue

        leftmost_x = valid_df[log2fc_col].min()
        rightmost_x = valid_df[log2fc_col].max()
        topmost_y = valid_df['-log10(p-value)'].max()

        for _, row in valid_df.iterrows():
            x = row[log2fc_col]
            y = row['-log10(p-value)']
            compound_id = row['Compounds ID']
            distance = np.sqrt((x - (leftmost_x if x < 0 else rightmost_x)) ** 2 + (y - topmost_y) ** 2)
            compound_distances[compound_id].append(distance)

    final_distances = {}
    for compound_id, distances in compound_distances.items():
        distances = np.array(distances)

        if len(distances) > 1 and np.var(distances) == 0:
            print(f"[DEBUG] Zero variance for compound {compound_id}, using simple mean.")
            final_distance = np.mean(distances)
        elif len(distances) > 1:
            final_distance = np.average(distances, weights=np.full_like(distances, 1 / np.var(distances)))
        else:
            final_distance = distances[0]

        final_distances[compound_id] = final_distance

    compound_distance_df = pd.DataFrame.from_dict(final_distances, orient='index', columns=['Total Distance'])

    mapping_df = pd.read_csv(file, encoding="ISO-8859-1", engine="python", sep=",", quotechar='"', on_bad_lines="skip")
    mapping_df.columns = clean_column_names(mapping_df.columns)

    # Clean and prepare mapping DataFrame
    required_cols = {'Compounds ID', 'Name', 'Formula', 'Calc. MW'}
    missing = required_cols - set(mapping_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in mapping file: {missing}")

    mapping = mapping_df[['Compounds ID', 'Name', 'Formula', 'Calc. MW']].drop_duplicates()

    # Sanitize quote pollution and fallback
    for col in ['Name', 'Formula', 'Calc. MW']:
        mapping[col] = mapping[col].apply(strip_nested_quotes)

    mapping['Name'] = mapping['Name'].fillna(mapping['Formula']).fillna(mapping['Calc. MW'])

    missing_in_mapping = set(compound_distance_df.index) - set(mapping['Compounds ID'])
    if missing_in_mapping:
        print(f"\n[INFO] {len(missing_in_mapping)} compounds in distance results not found in mapping.")

    compound_distance_named_df = compound_distance_df.merge(mapping, left_index=True, right_on='Compounds ID')
    # Sort and reorder final output columns
    final_result = compound_distance_named_df[['Compounds ID', 'Calc. MW', 'Name', 'Total Distance']]
    final_result_filtered = final_result[final_result['Total Distance'] < distance_threshold]

    if final_result_filtered.empty:
        print("\n[INFO] No compounds passed the distance threshold filter. Exiting early.")
        return pd.DataFrame(columns=['Compounds ID', 'Calc. MW', 'Name', 'Total Distance'])

    if len(final_result_filtered) > max_results:
        final_result_filtered = final_result_filtered.head(max_results)


    # Final validation: no names should be missing or malformed
    bad_names = final_result_filtered['Name'].apply(lambda x: pd.isna(x) or str(x).strip() in ('', '""', '``'))
    if bad_names.any():
        print(f"\n[ERROR] Found {bad_names.sum()} entries with unresolved names after all fallbacks.")
        print(final_result_filtered[bad_names])

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_colwidth', None)
    print("\nFinal filtered result:")
    print(final_result_filtered)

    final_result_filtered.to_csv("by_distance_named.csv", index=False)
    return final_result_filtered