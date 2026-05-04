# QC pipeline for a raw peptide x sample matrix:
# 1) Filter samples first: keep samples with >= 20 detected peptides
# 2) Filter peptides next: keep peptides detected in >= 20% of remaining samples
# 3) Rename sample columns: '1668_finished_1_1_WM_...' -> '1668_18_WM'
# 4) Save final matrix + summaries + removed lists

import re
# We use regex to pull structured pieces (TMA, row, col, region) out of long sample names.

import math
# We use math for ceiling so the peptide detection cutoff rounds up correctly.

import pandas as pd
# We use pandas to read the matrix from CSV and do table-style filtering and summaries.

from pathlib import Path
# We use Path to build file paths safely across folders.


# ===== 0. Helpers: basic safety checks =====
def check_labels_are_strings(df: pd.DataFrame) -> None:
    # This function checks that every row label and column label is a real, non-empty string, so downstream steps do not silently break.

    index_errors = []

    for title in df.index:
        # For each peptide ID in the index, confirm it is a non-empty string; otherwise store it as an error.
        if not isinstance(title, str) or title.strip() == "":
            index_errors.append(title)

    col_errors = []

    for title in df.columns:
        # For each sample name in the columns, confirm it is a non-empty string; otherwise store it as an error.
        if not isinstance(title, str) or title.strip() == "":
            col_errors.append(title)

    msgs = []

    if index_errors:
        # If any row labels are invalid, record a preview so the error message is actionable.
        msgs.append(
            "Row label errors (index, should be non-empty strings). "
            f"Examples: {index_errors[:20]} (total {len(index_errors)})"
        )

    if col_errors:
        # If any column labels are invalid, record a preview so the error message is actionable.
        msgs.append(
            "Column label errors (samples, should be non-empty strings). "
            f"Examples: {col_errors[:20]} (total {len(col_errors)})"
        )

    if msgs:
        # If any label issue exists, stop immediately so we do not generate corrupted QC output.
        raise ValueError("Label validation failed:\n" + "\n".join(msgs))


def coerce_matrix_numeric(df: pd.DataFrame) -> pd.DataFrame:
    # This function forces the whole matrix to numeric; if any entry is missing or cannot convert, it reports the exact peptide and sample and stops.

    out = df.copy()

    error_messages = []

    for col in out.columns:
        # For each sample column, attempt numeric conversion and detect any bad entries.
        original = out[col]
        converted = pd.to_numeric(original, errors="coerce")

        # If conversion creates NaN from something that was not NaN before, that means the value was non-numeric text or invalid.
        bad_convert_mask = converted.isna() & original.notna()

        # If the original already contains NaN, we treat that as invalid input for this QC pipeline.
        bad_missing_mask = original.isna()

        bad_mask = bad_convert_mask | bad_missing_mask

        if bad_mask.any():
            # If any cell in this column is invalid, collect peptide IDs and store a precise location + value for debugging.
            bad_peptides = out.index[bad_mask]
            for pep in bad_peptides:
                bad_val = out.at[pep, col]
                error_messages.append(
                    f"Bad value at peptide '{pep}', sample '{col}': {repr(bad_val)}"
                )

        out[col] = converted

    if error_messages:
        # If any invalid values exist anywhere in the matrix, stop and show a preview so you can fix the source file.
        preview = error_messages[:50]
        msg = "NaN or non-numeric entries found in matrix. Fix these in the source CSV/matrix:\n" + "\n".join(preview)
        if len(error_messages) > 50:
            msg += f"\n... and {len(error_messages) - 50} more."
        raise ValueError(msg)

    return out


def make_unique(names):
    # This function takes a list of names and ensures uniqueness by appending _2, _3, etc. when duplicates occur after renaming.

    seen = {}
    new = []

    for n in names:
        # For each proposed name, either keep it if new, or add a numeric suffix if we already saw it.
        if n not in seen:
            seen[n] = 1
            new.append(n)
        else:
            seen[n] += 1
            new.append(f"{n}_{seen[n]}")

    return new


# ===== 1. Helpers: renaming sample columns =====
def coord_to_core(row_idx, col_idx, n_cols=3, total_cores=18):
    # This function converts (row, col) grid coordinates into a core number using a fixed plate layout and a reversed numbering scheme.

    k = (row_idx - 1) * n_cols + col_idx
    core = total_cores + 1 - k

    return core


def relabel_column_to_short(col_name: str) -> str:
    # This function parses a long sample name, extracts TMA/row/col/region, converts (row,col) to a core number, then returns a short standardized label.

    m = re.match(
        r"^(?P<tma>\d+)_finished_(?P<row>\d+)_(?P<col>\d+)_(?P<region>WM|GM|LC).*?$",
        col_name
    )

    if not m:
        # If the name does not match the expected pattern, we do not rename it so we do not accidentally corrupt unknown formats.
        return col_name

    tma = m.group("tma")
    row = int(m.group("row"))
    col = int(m.group("col"))
    region = m.group("region")

    core = coord_to_core(row, col)

    return f"{tma}_{core}_{region}"


# ===== 2. Paths =====
data_dir = Path(r"C:\Users\tanji\Desktop\GM 40 Prec Thershold")
# This sets the project base folder so all outputs go to a consistent location.

summary_dir = data_dir / "Summary Folder"
summary_dir.mkdir(exist_ok=True)
# If the summary folder does not exist, create it so saving files will not fail.

result_dir = summary_dir
# This is the folder where all output summaries and the final matrix are written.

matrix_file = summary_dir / "GM_Raw_matrix.csv"
# This points to the input matrix file that will be QC-filtered.

if not matrix_file.exists():
    # If the input matrix is missing, stop here because nothing downstream can run correctly.
    raise FileNotFoundError(f"Raw matrix file not found at: {matrix_file}")

print(f"Loading raw matrix from: {matrix_file}")


# ===== 3. Load matrix =====
matrix_raw = pd.read_csv(matrix_file, index_col=0, encoding="latin1")
# This loads the matrix with peptide IDs as rows and sample names as columns.

matrix_raw.index = matrix_raw.index.astype(str).str.strip()
# We clean peptide IDs into consistent strings.

matrix_raw.columns = matrix_raw.columns.astype(str).str.strip()
# We clean sample names into consistent strings.

print(f"Raw matrix shape (peptides x samples): {matrix_raw.shape}")

check_labels_are_strings(matrix_raw)
# We validate that the matrix labels are usable for indexing and saving.

matrix_raw = coerce_matrix_numeric(matrix_raw)
# We ensure every intensity is numeric and fail fast if any bad entries exist.


# ===== 4. Presence definition =====
presence_threshold = 0.0
# This is the detection cutoff: values greater than this are considered "present".

presence = matrix_raw > presence_threshold
# This converts the intensity matrix into a boolean detected/not-detected matrix.


# ===== 5. Sample QC first =====
sample_peptide_hits = presence.sum(axis=0)
# For each sample, count how many peptides are detected.

sample_summary = pd.DataFrame({
    "sample": sample_peptide_hits.index,
    "peptide_hits": sample_peptide_hits.values
}).sort_values("peptide_hits", ascending=False)
# This builds a ranked table so you can see which samples have the most or fewest detected peptides.

sample_summary_file = result_dir / "raw_summary_samples_peptide_hits.xlsx"
sample_summary.to_excel(sample_summary_file, index=False)
# Save the per-sample detection summary for QC reporting.

min_sample_hits = 20
# Samples below this detection count will be removed.

good_samples = sample_peptide_hits.index[sample_peptide_hits >= min_sample_hits]
bad_samples = sample_peptide_hits.index[sample_peptide_hits < min_sample_hits]
# This splits samples into pass/fail groups based on the detection cutoff.

bad_samples_file = result_dir / "QC_removed_samples.txt"
with open(bad_samples_file, "w", encoding="utf-8") as f:
    # For all samples that failed the cutoff, write their original names to a text file for traceability.
    for s in bad_samples:
        f.write(str(s) + "\n")

matrix_samp = matrix_raw.loc[:, matrix_raw.columns.isin(good_samples)]
# This keeps only samples that pass the minimum detected-peptides requirement.

print(f"After sample filter shape: {matrix_samp.shape}")

if matrix_samp.shape[1] == 0:
    # If no samples remain after QC, stop because peptide filtering and saving would be meaningless.
    raise ValueError("All samples removed.")


# ===== 6. Peptide QC second =====
presence_samp = matrix_samp > presence_threshold
# Recompute detection after sample filtering because the remaining sample set changed.

n_samples_after = matrix_samp.shape[1]
# This is the number of samples that remain, which determines the peptide detection requirement.

min_fraction = 0.40
min_peptide_samples = int(math.ceil(min_fraction * n_samples_after))
# Convert the fraction cutoff into a required sample count, rounding up so the rule is strict.

peptide_sample_hits = presence_samp.sum(axis=1)
# For each peptide, count in how many remaining samples it is detected.

peptide_summary = pd.DataFrame({
    "Peptide_ID": peptide_sample_hits.index,
    "n_samples_positive": peptide_sample_hits.values
}).sort_values("n_samples_positive", ascending=False)
# This builds a ranked table so you can see which peptides are broadly detected versus rare.

peptide_summary_file = result_dir / "raw_summary_peptides_sample_hits.xlsx"
peptide_summary.to_excel(peptide_summary_file, index=False)
# Save the per-peptide detection summary for QC reporting.

good_peptides_mask = peptide_sample_hits >= min_peptide_samples
bad_peptides = peptide_sample_hits.index[~good_peptides_mask]
# This flags peptides that pass the detection frequency cutoff and records the ones to remove.

bad_peptides_file = result_dir / "QC_removed_peptides.txt"
with open(bad_peptides_file, "w", encoding="utf-8") as f:
    # For all peptides that failed the cutoff, write their IDs to a text file for traceability.
    for p in bad_peptides:
        f.write(str(p) + "\n")

matrix_qc = matrix_samp.loc[good_peptides_mask]
# This keeps only peptides that are detected in enough of the remaining samples.

print(f"After peptide filter shape: {matrix_qc.shape}")

if matrix_qc.shape[0] == 0:
    # If no peptides remain after QC, stop because the final matrix would be empty.
    raise ValueError("All peptides removed.")


# ===== 7. Rename columns =====
new_cols = [relabel_column_to_short(c) for c in matrix_qc.columns]
# For each original sample name, generate a short standardized label based on regex parsing.

new_cols = make_unique(new_cols)
# If renaming caused duplicates, append suffixes so column labels stay unique.

matrix_qc.columns = new_cols
# Apply the final renamed sample labels to the QC matrix.


# ===== 8. Final validation =====
check_labels_are_strings(matrix_qc)
# Confirm the final matrix labels are still valid after filtering and renaming.

matrix_qc = coerce_matrix_numeric(matrix_qc)
# Confirm filtering did not introduce any invalid numeric values.


# ===== 9. Save final matrix =====
final_matrix_file = result_dir / "FINAL_GM_QC_matrix.csv"
matrix_qc.to_csv(final_matrix_file)
# Save the final QC-filtered matrix so it is ready for downstream analysis.

print("Done!")
print(f"Final matrix shape: {matrix_qc.shape}")
