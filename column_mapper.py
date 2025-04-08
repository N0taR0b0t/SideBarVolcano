# column_mapper.py

import pandas as pd
import re
import json

def clean_key(s: str) -> str:
    s = re.sub(r'["\\]', '', s)
    s = re.sub(r'\s+', ' ', s)
    return s.strip()

def deep_clean(obj):
    if isinstance(obj, dict):
        return {k: deep_clean(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [deep_clean(i) for i in obj]
    elif isinstance(obj, str):
        v = obj.strip()
        v = v.replace('\\"', '')
        v = v.replace('"', '')
        v = v.replace('\\', '')
        v = re.sub(r'\s+', ' ', v)
        return v.strip()
    return obj

def extract_column_mapping(csv_file: str) -> dict:
    df = pd.read_csv(csv_file, encoding='latin1', nrows=0)
    columns = df.columns.tolist()

    fold_change_pattern = r"Log2 Fold Change: (.+)"
    p_value_pattern = r"(?<!Adj\.\s)P-value: (.+)"

    fold_change_cols = {
        clean_key(re.search(fold_change_pattern, col).group(1)): col.strip()
        for col in columns if re.search(fold_change_pattern, col)
    }
    p_value_cols = {
        clean_key(re.search(p_value_pattern, col).group(1)): col.strip()
        for col in columns if re.search(p_value_pattern, col)
    }

    matched_keys = set(fold_change_cols.keys()) & set(p_value_cols.keys())

    entries = []
    for key in matched_keys:
        entries.append({
            "fold_change_col": fold_change_cols[key],
            "p_value_col": p_value_cols[key],
            "title": f"Volcano Plot: {key}"
        })

    return {csv_file: deep_clean(entries)}

def save_mapping(mapping: dict, output_path: str):
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(mapping, f, indent=4, ensure_ascii=False)