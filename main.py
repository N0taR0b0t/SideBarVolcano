import os
import csv
import sys
from column_mapper import extract_column_mapping, save_mapping
from distCalc import compute_compound_distances
# from ModularVolcanos import generate_volcano_plot

def is_valid_csv(filename):
    return (
        filename.endswith(".csv") and
        filename.startswith("F_Re") and
        os.path.splitext(filename)[0].count("_") <= 1
    )

def is_truly_csv(filepath):
    """Check if the file is actually comma-delimited."""
    try:
        with open(filepath, 'r', newline='', encoding='latin1') as file:
            sample = file.read(4096)
            dialect = csv.Sniffer().sniff(sample, delimiters=[',', '\t'])
            return dialect.delimiter == ','
    except Exception as e:
        print(f"⚠️ Could not determine delimiter for file '{filepath}': {e}")
        return False

def main():
    csv_files = [f for f in os.listdir(".") if is_valid_csv(f)]

    # Check if all files are truly comma-delimited
    for csv_file in csv_files:
        if not is_truly_csv(csv_file):
            print(f"WARNING: File '{csv_file}' is NOT comma-delimited. It may be a TSV!")
            print("Please correct the delimiter or file format before proceeding.")
            sys.exit(1)

    if not csv_files:
        print("❌ No valid CSV files found in the current directory.")
        return

    for csv_file in csv_files:
        base_name = os.path.splitext(csv_file)[0]
        mapping_file = f"{base_name}_column_mapping.json"
        distance_file = f"{base_name}_by_distance_named.csv"
        #plot_file = f"{base_name}_volcano_plot.html"

        # Step 1: Generate mapping
        mapping = extract_column_mapping(csv_file)
        save_mapping(mapping, mapping_file)
        print(f"✅ Cleaned and saved mapping for '{csv_file}' to: {mapping_file}")

        # Step 2: Compute distances
        compute_compound_distances(csv_file, mapping_file, output_csv=distance_file)
        print(f"✅ Distance computation complete for '{csv_file}'.")

        # Step 3: Generate volcano plot
        #generate_volcano_plot(csv_file, mapping_file, distance_file, plot_file)
        #print(f"✅ Volcano plot saved to '{plot_file}'.")

if __name__ == "__main__":
    main()