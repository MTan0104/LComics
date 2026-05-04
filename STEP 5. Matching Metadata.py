import os
import re
import pandas as pd
import numpy as np


# ===== 0. Paths =====
os.chdir(r"C:\Users\tanji\Desktop\GM 40 Prec Thershold\Summary Folder")
# We set the working directory so reads and writes land in the same folder.


# ===== 1. Read matrix and metadata =====
matrix_path = "FINAL_GM_QC_matrix.csv"
meta_path = "FinalDatabase.csv"

matrix_df = pd.read_csv(matrix_path, dtype=str)
# We read as strings first so parsing column names is safe.

meta_df = pd.read_csv(meta_path)
# We read metadata normally since it contains mixed types.

if "Peptide_ID" not in matrix_df.columns:
    raise ValueError("Matrix must contain a column named 'Peptide_ID'.")
# If Peptide_ID is missing, we cannot separate sample columns from peptide IDs cleanly.

required_meta_cols = ["TMA", "Core", "TissueDiagnosisID"]
missing_meta_cols = [c for c in required_meta_cols if c not in meta_df.columns]
if missing_meta_cols:
    raise ValueError(f"Metadata is missing required columns: {missing_meta_cols}")
# The join logic depends on TMA/Core and TissueDiagnosisID for TissueRegionID.

print(f"Loaded matrix: {matrix_path}  shape={matrix_df.shape}")
print(f"Loaded metadata: {meta_path}  shape={meta_df.shape}")


# ===== 2. Clean and validate metadata keys =====
meta_df["TMA"] = pd.to_numeric(meta_df["TMA"], errors="raise").astype(int)
meta_df["Core"] = pd.to_numeric(meta_df["Core"], errors="raise").astype(int)
# We force TMA/Core to int so matches are exact and not string-vs-int mismatches.

meta_df["TissueDiagnosisID"] = meta_df["TissueDiagnosisID"].astype(str).str.strip()
# We make TissueDiagnosisID a clean string because we will concatenate it with Region.

if meta_df[["TMA", "Core"]].isna().any().any():
    raise ValueError("Metadata has missing values in TMA or Core. Fix metadata first.")
# Missing keys make matching ambiguous and should not proceed.

dup_core_keys = meta_df.duplicated(subset=["TMA", "Core"], keep=False)
if dup_core_keys.any():
    bad = meta_df.loc[dup_core_keys, ["TMA", "Core", "TissueDiagnosisID"]].sort_values(["TMA", "Core"])
    raise ValueError(
        "Metadata has duplicate rows for the same (TMA, Core). Matching would be ambiguous.\n"
        "Examples:\n" + bad.head(30).to_string(index=False)
    )
# If (TMA, Core) repeats, one matrix column could map to multiple metadata rows, which is not allowed.


# ===== 3. Parse matrix sample columns into (TMA, Core, Region) =====
sample_cols = [c for c in matrix_df.columns if c != "Peptide_ID"]
# These are the columns that should match metadata after thresholding.

if len(sample_cols) == 0:
    raise ValueError("Matrix has no sample columns besides Peptide_ID.")
# If no sample columns exist, there is nothing to align.

pattern = re.compile(r"^(?P<tma>\d+)_(?P<core>\d+)_(?P<region>GM|WM|LC)$")
# This defines the required naming format for sample columns.

parsed = []
bad_cols = []

for c in sample_cols:
    # For each sample column name, extract TMA/Core/Region from the string.
    m = pattern.match(str(c).strip())
    if not m:
        bad_cols.append(c)
        continue
    parsed.append({
        "ColumnName": c,
        "TMA": int(m.group("tma")),
        "Core": int(m.group("core")),
        "Region": m.group("region"),
    })

if bad_cols:
    # If any columns do not match the expected format, stop immediately so you do not silently mislabel samples.
    preview = bad_cols[:30]
    raise ValueError(
        "Some matrix columns do not match expected format 'TMA_Core_Region' (example: 1668_1_GM).\n"
        f"Examples: {preview}\n"
        f"Total bad columns: {len(bad_cols)}"
    )

region_info = pd.DataFrame(parsed)
# This table represents the actual samples left in the matrix.

regions_in_matrix = sorted(region_info["Region"].unique().tolist())
print(f"Regions present in matrix: {regions_in_matrix}")
# This makes the script automatically work for GM-only, WM-only, LC-only, or all-region matrices.


# ===== 4. Triplicate metadata into region-level rows =====
meta_base = meta_df.copy()
# One row per (TMA, Core) in your original metadata.

meta_long = meta_base.loc[meta_base.index.repeat(len(regions_in_matrix))].copy()
# We repeat each metadata row once per region, but only for regions that actually exist in the matrix.

meta_long["Region"] = np.tile(regions_in_matrix, len(meta_base))
# We assign the Region labels in the same repeated pattern.

meta_long["TissueRegionID"] = meta_long["TissueDiagnosisID"].astype(str) + "_" + meta_long["Region"]
# This creates your desired label like 93760_GM, 93760_WM, 93760_LC.


# ===== 5. Match matrix columns to metadata rows =====
meta_for_matrix = region_info.merge(
    meta_long,
    how="left",
    on=["TMA", "Core", "Region"],
)
# For each matrix column, find the corresponding metadata row by matching (TMA, Core, Region).

if meta_for_matrix["TissueDiagnosisID"].isna().any():
    # If any matrix sample cannot find metadata, stop and show which ones failed so you can diagnose core numbering or missing rows.
    missing = meta_for_matrix.loc[
        meta_for_matrix["TissueDiagnosisID"].isna(),
        ["ColumnName", "TMA", "Core", "Region"]
    ]
    raise ValueError(
        "Some matrix columns did not find matching metadata using (TMA, Core, Region).\n"
        "This usually means metadata is missing that (TMA, Core), or your Core numbering convention differs between files.\n"
        "Examples:\n" + missing.head(30).to_string(index=False)
    )

# Move TissueRegionID to column 2 for readability.
cols = meta_for_matrix.columns.tolist()
if "TissueRegionID" in cols:
    cols.remove("TissueRegionID")
    cols.insert(1, "TissueRegionID")
    meta_for_matrix = meta_for_matrix[cols]

meta_for_matrix.to_csv("metadata_MATCHED.csv", index=False)
# This file contains metadata aligned to matrix columns (one row per matrix sample column).

print("Saved: metadata_MATCHED.csv")


# ===== 6. Find metadata rows that are missing from the matrix =====
all_keys = meta_long[["TMA", "Core", "Region"]].drop_duplicates()
# This represents every (TMA, Core, Region) that could exist based on metadata expansion.

used_keys = region_info[["TMA", "Core", "Region"]].drop_duplicates()
# This represents the (TMA, Core, Region) that actually exist in the matrix columns.

missing_keys = all_keys.merge(
    used_keys,
    on=["TMA", "Core", "Region"],
    how="left",
    indicator=True,
)
missing_keys = missing_keys[missing_keys["_merge"] == "left_only"].drop(columns="_merge")
# These are region-expanded metadata rows that do not appear in the matrix (likely thresholded out).

missing_meta = meta_long.merge(
    missing_keys,
    on=["TMA", "Core", "Region"],
    how="inner",
)
# This returns the full metadata for all missing region-level rows.

missing_meta.to_csv("metadata_removed_rows.csv", index=False)
# This is your record of what metadata rows were not used because the matrix does not contain those samples.

print("Saved: metadata_removed_rows.csv")
print("Done! Your New Metadata file is saved.")
