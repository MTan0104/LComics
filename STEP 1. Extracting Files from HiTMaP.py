from pathlib import Path
import shutil

# ===== 1. Root folder where all the "XXXX finished" folders live =====
root = Path(r"E:\Finished Core") # This is the top-level folder. Inside this folder are many subfolders like "1173 finished", "1174 finished",etc.


# ===== 2. Destination folder on Desktop =====
desktop = Path.home() / "Desktop" # Path.home() finds your user home folder. Then we attach "Desktop" to it. This makes it portable across machines.
dest = desktop / "All_Peptide_Summaries" # Inside Desktop, create (or use) a folder called "All_Peptide_Summaries".
dest.mkdir(parents=True, exist_ok=True) # If the folder does not exist, create it. Parents=True allows creation of parent folders if needed. Exist_ok=True prevents crashing if it already exists.
print(f"Saving files to: {dest}")


# ===== 3. Find and copy all Peptide_Summary Excel files =====
for peptide_file in root.rglob("Peptide_Summary*.csv"):

    summary_folder = peptide_file.parent
    sample_folder = summary_folder.parent
    core_folder = sample_folder.parent

    core_name = core_folder.name.replace(" ", "_")
    sample_name = sample_folder.name.replace(" ", "_")

    new_name = f"{core_name}_{sample_name}_Peptide_Summary.csv"
    target = dest / new_name

    print(f"Copying:\n {peptide_file}\n -> {target}")

    shutil.copy2(peptide_file, target)

print("Done!")
