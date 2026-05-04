from pathlib import Path
import pandas as pd

# ====== SETTINGS ======
in_dir = Path(r"C:\Users\tanji\Desktop\MS Analysis 2026\HiTMaP File Sorting\All_Protein_Files")
out_csv = in_dir / "ALL_UNIQUE_FEATURES_LIST.csv"
conflict_txt = in_dir / "MZ_PEPTIDE_CONFLICT_REPORT.txt"
skipped_txt = in_dir / "SKIPPED_FILES_REPORT.txt"

pattern = "*.csv"

protein_col = "Protein"
mz_col = "mz"
peptide_col = "Peptide"
mz_align_col = "mz_align"
desc_col = "desc"
mod_col = "Modification"
adduct_col = "adduct"

needed_cols = [protein_col, mz_col, peptide_col, mz_align_col, desc_col, mod_col, adduct_col]
# ======================

all_rows = []
skipped = []

files = sorted(in_dir.glob(pattern))

for f in files:
    # Skip outputs/logs we generate in the same folder
    if f.name in {out_csv.name, conflict_txt.name, skipped_txt.name, "copy_log.csv"}:
        continue

    try:
        df = pd.read_csv(f, dtype=str, low_memory=False)
    except UnicodeDecodeError:
        df = pd.read_csv(f, dtype=str, low_memory=False, encoding="latin1")
    except Exception as e:
        skipped.append((f.name, f"read_failed: {e}"))
        continue

    missing_cols = [col for col in needed_cols if col not in df.columns]
    if missing_cols:
        skipped.append((f.name, f"missing_cols: {missing_cols}"))
        continue

    tmp = df[needed_cols].copy()
    tmp["source_file"] = f.name

    # Clean spaces
    for col in needed_cols:
        tmp[col] = tmp[col].astype(str).str.strip()

    # Remove empty / nan-like mz rows
    tmp = tmp[(tmp[mz_col] != "") & (tmp[mz_col].str.lower() != "nan")]

    # Clean other columns
    for col in [protein_col, peptide_col, mz_align_col, desc_col, mod_col, adduct_col]:
        tmp[col] = tmp[col].replace({"nan": "", "None": "", "<NA>": ""})

    # Remove exact duplicate rows within this file
    tmp = tmp.drop_duplicates()

    all_rows.append(tmp)

if not all_rows:
    raise RuntimeError("No valid files found with required columns. Check folder path and file pattern.")

master = pd.concat(all_rows, ignore_index=True)

# --------------------------------------------------
# Step 1: report mz values linked to multiple peptides
# --------------------------------------------------
conflict_lines = []

for mz_value, grp in master.groupby(mz_col, dropna=False):
    peptide_values = sorted(set(x.strip() for x in grp[peptide_col].fillna("") if x.strip() != ""))

    if len(peptide_values) > 1:
        conflict_lines.append("=" * 100)
        conflict_lines.append(f"m/z CONFLICT: {mz_value}")
        conflict_lines.append(f"Distinct peptides found: {len(peptide_values)}")
        conflict_lines.append("Peptides:")
        for pep in peptide_values:
            conflict_lines.append(f"  - {pep}")
        conflict_lines.append("")
        conflict_lines.append("Rows involved:")
        for _, row in grp[[mz_col, protein_col, peptide_col, mod_col, adduct_col, mz_align_col, desc_col, "source_file"]].drop_duplicates().iterrows():
            conflict_lines.append(
                f"  file={row['source_file']} | Protein={row[protein_col]} | "
                f"Peptide={row[peptide_col]} | Mod={row[mod_col]} | Adduct={row[adduct_col]} | "
                f"mz_align={row[mz_align_col]} | desc={row[desc_col]}"
            )
        conflict_lines.append("")

if conflict_lines:
    with open(conflict_txt, "w", encoding="utf-8") as fh:
        fh.write("\n".join(conflict_lines))
    print(f"Conflict report written: {conflict_txt}")
else:
    with open(conflict_txt, "w", encoding="utf-8") as fh:
        fh.write("No m/z to peptide conflicts found.\n")
    print(f"No conflicts found. Report written: {conflict_txt}")

# --------------------------------------------------
# Step 2: keep all unique feature rows
# --------------------------------------------------
# Keep unique combinations including modification and adduct
final = master.drop_duplicates(
    subset=[protein_col, mz_col, peptide_col, mod_col, adduct_col, mz_align_col, desc_col],
    keep="first"
).copy()

# Sort numerically by mz if possible
final["_mz_num"] = pd.to_numeric(final[mz_col], errors="coerce")
final = final.sort_values(
    by=["_mz_num", mz_col, peptide_col, protein_col],
    na_position="last"
).drop(columns=["_mz_num", "source_file"])

# Save final output
final.to_csv(out_csv, index=False)

print(f"Done. Unique feature rows: {len(final)}")
print(f"Saved feature table: {out_csv}")

# --------------------------------------------------
# Step 3: write skipped file report
# --------------------------------------------------
if skipped:
    with open(skipped_txt, "w", encoding="utf-8") as fh:
        for name, reason in skipped:
            fh.write(f"{name}\t{reason}\n")
    print(f"Some files were skipped. See: {skipped_txt}")