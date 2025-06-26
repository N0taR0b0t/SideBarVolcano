import os
import csv
import sys
import re
import pandas as pd
from column_mapper import extract_column_mapping, save_mapping
from distCalc import compute_compound_distances

# ─────────────────────────────────────────────────────────────
IDENTIFIER = "SysPlas"
REQUIRED_COLUMN = "Compounds ID"
ENCODING = "latin1"
DELIMITER = ","

# ─────────────────────────────────────────────────────────────
def is_valid_csv(filename):
    filename = filename
    return (
        filename.endswith(".csv")
        and filename.startswith(IDENTIFIER)
        and os.path.splitext(filename)[0].count("_") <= 1
    )

def sanitize_quotes_and_commas(raw_line: str, expected_fields: int) -> str or None:
    try:
        fields = list(csv.reader([raw_line], delimiter=DELIMITER))[0]
    except Exception:
        return None

    if len(fields) == expected_fields:
        return raw_line

    cleaned = []
    current = ""
    in_quotes = False
    for part in raw_line.strip().split(DELIMITER):
        part = part.strip()
        if part.startswith('"') and not part.endswith('"'):
            in_quotes = True
            current = part
        elif in_quotes:
            current += f",{part}"
            if part.endswith('"'):
                cleaned.append(current)
                in_quotes = False
                current = ""
        else:
            cleaned.append(part if not current else current + f",{part}")
            current = ""

    cleaned = [f.strip().strip('"') for f in cleaned]
    if len(cleaned) < expected_fields:
        return None
    elif len(cleaned) > expected_fields:
        cleaned = cleaned[:expected_fields]

    return DELIMITER.join([f'"{f}"' for f in cleaned])

def clean_csv_file(input_path: str, output_path: str) -> int:
    with open(input_path, "r", encoding=ENCODING) as infile:
        raw_lines = infile.readlines()

    header = raw_lines[0].strip()
    expected_field_count = len(list(csv.reader([header], delimiter=DELIMITER))[0])
    cleaned_lines = [header + "\n"]
    num_fixed = 0
    num_skipped = 0

    for lineno, line in enumerate(raw_lines[1:], start=2):
        line = line.strip()
        cleaned_line = sanitize_quotes_and_commas(line, expected_field_count)
        if cleaned_line is None:
            print(f"⚠️ Line {lineno} could not be fixed. Skipping.")
            num_skipped += 1
            continue
        field_count = len(list(csv.reader([cleaned_line], delimiter=DELIMITER))[0])
        if field_count != expected_field_count:
            print(f"⚠️ Line {lineno} still malformed after fix ({field_count} fields). Skipping.")
            num_skipped += 1
            continue
        cleaned_lines.append(cleaned_line + "\n")
        num_fixed += 1

    with open(output_path, "w", encoding=ENCODING, newline="") as outfile:
        outfile.writelines(cleaned_lines)

    print(f"ℹ️ Skipped {num_skipped} malformed rows.")
    return num_fixed

def validate_compound_ids(df: pd.DataFrame, source_file: str):
    if REQUIRED_COLUMN not in df.columns:
        raise ValueError(f"❌ Missing required column '{REQUIRED_COLUMN}' in file: {source_file}")

    df[REQUIRED_COLUMN] = df[REQUIRED_COLUMN].astype(str).str.strip()
    num_missing = df[REQUIRED_COLUMN].replace("", pd.NA).isna().sum()

    if num_missing > 0:
        raise ValueError(f"❌ Found {num_missing} missing '{REQUIRED_COLUMN}' values in '{source_file}'.")

def load_clean_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, encoding=ENCODING)
    df.columns = df.columns.str.strip().str.replace('"', "", regex=False)
    df = df.applymap(lambda x: str(x).strip().replace('"', "") if pd.notnull(x) else x)
    return df

# ─────────────────────────────────────────────────────────────
def main():
    files_in_dir = os.listdir(".")
    csv_files = [f for f in files_in_dir if is_valid_csv(f)]

    if not csv_files:
        print("❌ No valid CSV files found.")
        sys.exit(1)

    for csv_file in csv_files:
        print(f"\n🔍 Processing: {csv_file}")

        base_name = os.path.splitext(csv_file)[0]
        cleaned_file = base_name + ".csv"
        mapping_file = f"{base_name}_column_mapping.json"
        distance_file = f"{base_name}_by_distance_named.csv"

        try:
            rows_kept = clean_csv_file(csv_file, cleaned_file)
            print(f"✅ Cleaned {rows_kept} rows → {cleaned_file}")

            df = load_clean_csv(cleaned_file)
            validate_compound_ids(df, cleaned_file)
            print(f"✅ '{REQUIRED_COLUMN}' column is valid.")

            mapping = extract_column_mapping(cleaned_file)
            save_mapping(mapping, mapping_file)
            print(f"✅ Mapping saved → {mapping_file}")

            compute_compound_distances(cleaned_file, mapping_file, output_csv=distance_file)
            print(f"✅ Distance file written → {distance_file}")

        except Exception as e:
            print(f"\n🚨 Failed on file '{csv_file}': {e}")
            sys.exit(1)

# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
