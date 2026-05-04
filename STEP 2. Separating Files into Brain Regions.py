# Step 1: Separate Peptide_Summary CSV files into folders by brain region (GM, WM, LC)

import shutil  
# Hey Python, I need a toolbox that can copy and move files around.

from pathlib import Path  
# Hey Python, I want to work with folders and file paths in a clean way.
# Path is like a smart folder/file object instead of raw strings.

# ===== 1. Base folder with all Peptide_Summary files =====
data_dir = Path(r"C:\Users\tanji\Desktop\GM Proteomic Analysis 2026")
# Hey Python, this is the main folder where ALL my CSV files live.
# r"" means raw string so backslashes are treated normally.
# Path(...) turns this into a Path object that Python can navigate.

# ===== 2. Destination folders =====
dest_base = data_dir / "ByRegion"
# Hey Python, inside my main folder, create (or refer to) a subfolder called "ByRegion".
# The "/" operator joins paths safely.

dest_gm = dest_base / "GM"
# Inside ByRegion, this will be the folder for GM files.

dest_wm = dest_base / "WM"
# Inside ByRegion, this will be the folder for WM files.

dest_lc = dest_base / "LC"
# Inside ByRegion, this will be the folder for LC files.

for d in [dest_gm, dest_wm, dest_lc]:
    d.mkdir(parents=True, exist_ok=True)
# Hey Python, for each of these folders:
# If it does not exist, create it.
# parents=True means create parent folders if needed.
# exist_ok=True means do NOT crash if the folder already exists.

print(f"Source folder: {data_dir}")
# Print where we are reading files from.

print(f"Destination base folder: {dest_base}")
# Print where we are sending files to.

# ===== 3. Loop through all Peptide_Summary CSV files =====
pattern = "*Peptide_Summary*.csv"
# Hey Python, I only care about CSV files that contain "Peptide_Summary" in their name.
# The * means "anything before or after".

files = sorted(data_dir.glob(pattern))
# Hey Python, look inside data_dir.
# Find all files that match the pattern.
# Sort them alphabetically.
# Store them in a list called "files".

print(f"Found {len(files)} Peptide_Summary CSV files.")
# Count how many matching files we found and print it.

unmatched = []
# This will store files that do NOT match GM/WM/LC.
# Think of this as a warning list.

for f in files:
    # For each file in the list...

    name = f.name
    # Get just the filename (not the full path).

    upper_name = name.upper()
    # Convert the filename to uppercase.
    # This avoids problems like gm vs GM vs Gm.
    # Now matching is case-insensitive.

    if "_GM_" in upper_name:
        dest_folder = dest_gm
        # If filename contains "_GM_", send it to GM folder.

    elif "_WM_" in upper_name:
        dest_folder = dest_wm
        # If filename contains "_WM_", send it to WM folder.

    elif "_LC_" in upper_name:
        dest_folder = dest_lc
        # If filename contains "_LC_", send it to LC folder.

    else:
        # If none of those patterns matched...
        unmatched.append(name)
        # Add this filename to the unmatched list.

        print(f"  [SKIP] Could not detect region for: {name}")
        # Print a warning that we skipped it.

        continue
        # Skip the rest of the loop and move to the next file.

    target = dest_folder / name
    # Build the full destination path for the copy.
    # Example: ByRegion/GM/filename.csv

    print(f"Copying {name} -> {target}")
    # Tell us what is being copied and where.

    shutil.copy2(f, target)
    # Actually copy the file.
    # copy2 keeps metadata (like modification time).

print("\nDone separating files by region.")
# Finished processing all files.

print(f"GM files folder: {dest_gm}")
print(f"WM files folder: {dest_wm}")
print(f"LC files folder: {dest_lc}")
# Show final folder locations.

if unmatched:
    # If there are files that did not match...

    print("\nFiles that did NOT match GM/WM/LC pattern:")
    for n in unmatched:
        print("  -", n)
    # Print each unmatched filename.

else:
    print("\nAll files matched a region.")
    # If unmatched list is empty, everything was classified successfully.
