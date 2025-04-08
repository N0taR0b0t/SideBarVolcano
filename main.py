from column_mapper import extract_column_mapping, save_mapping
from distCalc import compute_compound_distances
from ModularVolcanos import generate_volcano_plot

csv_file = "ReSpleen.csv"
output_file = "column_mapping.json"

# Step 1: Generate mapping
mapping = extract_column_mapping(csv_file)
save_mapping(mapping, output_file)
print(f"✅ Cleaned and saved mapping for '{csv_file}' to: {output_file}")

# Step 2: Compute distances
compute_compound_distances(csv_file)
print("✅ Distance computation complete.")

# Step 3: Generate volcano plot
generate_volcano_plot()