# organizing all data into a matrix format from the raw csv files generated from HiTMaP

import pandas as pd
# Hey Python, I need pandas to read CSV files and reshape tables.

from pathlib import Path
# Hey Python, I want to handle folders and file paths safely using Path objects.

# ===== 1. Paths =====
data_dir = Path(r"C:\Users\tanji\Desktop\GM 40 Prec Thershold")
# Hey Python, this is the folder that contains all my GM Peptide_Summary CSV files.

output_file = data_dir / "GM_Raw_matrix.csv"
# Hey Python, this is the output file name I plan to use (for now) inside the data folder.

files = sorted(data_dir.glob("*.csv"))
# Hey Python, look inside data_dir and collect every .csv file.
# Sort them alphabetically so the run is reproducible.

print(f"Found {len(files)} files.")
# Print how many CSV files were found.

all_rows = []
# This will store small tables from each sample, all stacked together later (long format).

skipped_files = []   # <--- will store files that are missing columns
# This will store warnings about files that cannot be processed.

# ===== 2. Loop through every CSV file =====
for f in files:
    # Hey Python, process each CSV file one by one.

    sample_name = f.stem
    # Take the filename without extension as the sample name.
    # Example: "ABC123.csv" becomes "ABC123".

    print(f"Processing {sample_name} ...")
    # Print which sample is being processed.

    df = pd.read_csv(f)
    # Read the CSV into a pandas DataFrame.

    df.columns = df.columns.str.strip()
    # Clean column names by removing leading/trailing spaces.
    # This prevents errors if a column is named "Protein " instead of "Protein".

    required = {"Protein", "mz", "Peptide", "Modification", "adduct", "Intensity"}
    # Define the set of columns that MUST exist for this pipeline to work.

    missing = required - set(df.columns)
    # Compute which required columns are missing in this file.
    # This is set subtraction: required minus actual columns.

    if missing:
        # If any required columns are missing, we do not trust this file.

        msg = f"{f.name} (missing columns: {sorted(missing)})"
        # Build a readable message that lists which columns are missing.

        print("  -> SKIP:", msg)
        # Print a clear skip warning.

        skipped_files.append(msg)
        # Save this skip message so we can export a report later.

        continue
        # Stop processing this file and move to the next one.

    df["Modification"] = df["Modification"].fillna("")
    # Replace missing (NaN) modification values with an empty string.
    # This avoids "nan" showing up in IDs later.

    df["adduct"] = df["adduct"].fillna("")
    # Replace missing (NaN) adduct values with an empty string for the same reason.

    def make_id(row):
        # Hey Python, this function will build a unique string ID for each peptide feature.

        parts = [
            str(row["Protein"]).strip(),
            # Convert Protein to string and strip spaces.

            str(row["mz"]).strip(),
            # Convert mz to string and strip spaces.

            str(row["Peptide"]).strip(),
            # Convert peptide sequence to string and strip spaces.
        ]

        mod = str(row["Modification"]).strip()
        # Get modification text as a clean string.

        add = str(row["adduct"]).strip()
        # Get adduct text as a clean string.

        if mod:
            parts.append(mod)
            # Only include modification if it is not empty.

        parts.append(add)
        # Always include adduct, even if it is empty.
        # This forces consistent ID structure.

        return "_".join(parts)
        # Join all pieces with underscores to form the final Peptide_ID.

    df["Peptide_ID"] = df.apply(make_id, axis=1)
    # Create a new column by applying make_id to each row.
    # axis=1 means apply across rows (one row at a time).

    tmp = df.loc[:, ["Peptide_ID", "Intensity"]].copy()
    # Extract only the columns we care about for the matrix.
    # .copy() avoids pandas warnings and prevents accidental edits to df.

    tmp.columns = ["Peptide_ID", "intensity"]
    # Rename columns to a consistent naming scheme.

    tmp["sample"] = sample_name
    # Add a new column that labels which sample this intensity came from.

    all_rows.append(tmp)
    # Save this sample's long-format table into the list for later concatenation.

# ===== 5. Check =====
if not all_rows:
    # If nothing was collected, that means every file got skipped or there were no files.

    print("\nNo valid files were processed.")
    # Tell the user clearly.

    raise SystemExit
    # Stop the script immediately, because there is nothing to build a matrix from.

# ===== 6. Combine and pivot =====
long_df = pd.concat(all_rows, ignore_index=True)
# Stack all sample tables into one big long-format DataFrame.
# ignore_index=True resets row numbering from 0..N-1.

matrix = long_df.pivot_table(
    # Reshape from long format to wide matrix format.

    index="Peptide_ID",
    # Rows are peptide features.

    columns="sample",
    # Columns are sample names.

    values="intensity",
    # Each cell is the intensity value.

    fill_value=0,
    # If a peptide is missing in a sample, fill that entry with 0.
)

matrix.columns.name = None
# Remove the columns axis name ("sample") so the output CSV header is cleaner.

matrix = matrix.sort_index(axis=0).sort_index(axis=1)
# Sort rows by Peptide_ID and sort columns by sample name for consistent ordering.

# ===== Create "Summary Folder" and save there =====
summary_dir = data_dir / "Summary Folder"
# Define the folder where outputs and reports will be saved.

summary_dir.mkdir(exist_ok=True)
# Create the folder if it does not exist.
# exist_ok=True prevents errors if the folder already exists.

print(f"\nSummary folder created at: {summary_dir}")
# Print where the outputs will be saved.

output_file = summary_dir / "GM_Raw_matrix.csv"
# Update the output_file path so it saves into the Summary Folder.

matrix.to_csv(output_file)
# Save the final matrix as a CSV file.

print("\nDone!")
# Print completion message.

print(f"Saved matrix to: {output_file}")
# Print the output location.

# ===== Save skipped files list into Summary Folder =====
if skipped_files:
    # If any files were skipped, export a skip report.

    skipped_file_path = summary_dir / "GM_Raw_matrix_skipped_files.csv"
    # Define the path for the skip report.

    pd.DataFrame(skipped_files, columns=["skipped_file"]).to_csv(skipped_file_path, index=False)
    # Convert skip messages into a one-column table and save it as CSV.
    # index=False prevents pandas from writing an extra index column.

    print(f"\nSaved skipped file list to: {skipped_file_path}")
    # Print where the skip report was saved.

else:
    print("\nNo files were skipped.")
    # If the skip list is empty, everything passed the column checks.
