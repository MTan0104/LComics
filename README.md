# LComics
Code for data analysis of LC proteomic and glycomic data.

This project investigates region-specific molecular features in the human locus coeruleus (LC) using spatial proteomics and quantitative imaging within a human cohort of 126 individuals without clinically diagnosed neurodegenerative disease. Each tissue microarray (TMA) core was histologically annotated and manually segmented into gray matter (GM), white matter (WM), and LC neuron–enriched regions using QuPath, followed by spatial registration to MALDI mass spectrometry imaging (MALDI-MSI) data processed in SCiLS Lab.

Region-specific mass spectrometry data were exported as normalized spectral files (imzML and ibd) under total ion count normalization, enabling quantitative comparison of molecular signals across spatially defined tissue compartments. Peptide features were subsequently identified using the **HiTMaP** ([GitHub](https://github.com/MASHUOA/HiTMaP) | [Paper](https://doi.org/10.1038/s41467-021-23461-w)) computational framework, which maps MALDI-MSI spectra to in silico peptide libraries while controlling for false discovery.

**This repository implements a downstream computational workflow following HiTMaP output, including large-scale matrix reconstruction, quality filtering, and metadata integration**. The pipeline consolidates region-resolved peptide features into structured matrices, applies empirically defined thresholds to improve robustness, and aligns molecular data with sample-level metadata containing neuropathological and demographic variables such as Braak stage, Thal phase, and APOE genotype.

Together, this workflow provides a reproducible framework for transforming spatial proteomics outputs into analysis-ready datasets, enabling systematic investigation of region-specific molecular alterations associated with early pathological processes in the human brain.

## Step 1 – Extracting Files from HiTMaP
[Extracting Files from HiTMaP.py](https://github.com/user-attachments/files/27374387/Extracting.Files.from.HiTMaP.py)

This step collects all HiTMaP-generated Peptide_Summary.csv files from a nested directory structure and consolidates them into a single standardized folder. Each file is renamed during extraction to preserve its origin (core and sample identity) and to prevent filename collisions.

HiTMaP outputs are typically stored across multiple subdirectories corresponding to different cores and samples. This step performs a recursive search to locate all relevant peptide summary files and reorganizes them into a unified dataset, which serves as the input for downstream processing.

Input:
- Root directory containing multiple "XXXX finished" folders
  (each containing sample-level subfolders and summary outputs)

Process:
1. Recursively search all subfolders for Peptide_Summary.csv files
2. Extract core and sample identifiers from folder hierarchy
3. Rename files into a standardized format
4. Copy all files into a single destination folde

❗**Key assumptions (must be understood before use)**
1. Folder hierarchy structure
The code loop inside folder that been setted up by HitMaP:
```python
summary_folder = peptide_file.parent
sample_folder = summary_folder.parent
core_folder = sample_folder.parent

# Root/
# ├── XXXX finished/
# │   ├── Sample/
# │   │   └── Summary/
# │   │       └── Peptide_Summary.csv
```
The extraction logic depends on this structure, for example:
```python
#  E:\Finished Core #project specific naming
#  ├── 1173 finished #project specific naming
#  │   ├── Sample A
#  │   │   └── Summary
#  │   │       └── Peptide_Summary.csv
#  │   ├── Sample B
#  │       └── Summary
#  │           └── Peptide_Summary.csv
#  ├── 1174 finished
#      ├── Sample C
#          └── Summary
#              └── Peptide_Summary.csv
#  ......
```
For each file, You want to build:
```python
new_name = f"{core_name}_{sample_name}_Peptide_Summary.csv"
```
2. Folder names encode metadata
The script assumes:
-Core identity is stored in the core folder name
-Sample identity is stored in the sample folder name
```python
core_name = core_folder.name.replace(" ", "_")
sample_name = sample_folder.name.replace(" ", "_")
```
3. Output naming convention

Files are renamed as: [Core]_[Sample]_Peptide_Summary.csv. This naming is specific to this project, which should be adjust for other cases.

4. Destination location

```python
dest = Path.home() / "Desktop" / "All_Peptide_Summaries"
```
The output path will be the input path for **Step 2 - Separating Files into Brain Regions**

If a comprehensive feature set without region-specific separation is required (e.g., for global feature enumeration or cross-region comparison), please refer to [Step 6 (Optional) – Global Feature Extraction for Validation](#step-6--Global-Feature-Extraction-for-Validation), where all HiTMaP-derived features are consolidated into a unified dataset prior to region-specific processing.

## Step 2 - Separating Files into Brain Regions
[Separating Files into Brain Regions.py](https://github.com/user-attachments/files/27374390/Separating.Files.into.Brain.Regions.py)


This step organizes raw HiTMaP output files (Peptide_Summary.csv) into region-specific folders based on brain region annotations embedded in file names. All files are initially placed in a single directory and then automatically sorted into Gray Matter (GM), White Matter (WM), and Locus Coeruleus (LC) subfolders.

Input:
- Folder containing all HiTMaP Peptide_Summary.csv files from Step 1 Extracting Files from HiTMaP. 

Process:
- Detect region from filename
- Create region-specific folders (GM, WM, LC)
- Copy files into corresponding folders

Output:
ByRegion
- Gray Matter files
- White Matter files
- Locus Coeruleus files

❗**Key assumptions (must be understood before use)**
The script scans each file name for region identifiers and copies the files into corresponding directories, creating a structured dataset for downstream region-specific analysis.

```python
_GM_
_WM_
_LC_
```
Your file has to be named with those region identifiers or you can alter it for your project 

Before running the script, you must update the file path:
```python
data_dir = Path(r"YOUR_PATH_TO_PEPTIDE_SUMMARY_FILES")

#Example:

data_dir = Path(r"D:\LC_project\HiTMaP_outputs")
```

## Step 3 – Matrix Construction from HiTMaP Output
[Matrix Construction from HiTMaP Output.py](https://github.com/user-attachments/files/27371770/Matrix.Construction.from.HiTMaP.Output.py)

This step reconstructs a unified quantitative matrix from region-specific HiTMaP output files by transforming individual sample-level peptide summaries into a structured feature-by-sample matrix. Each input CSV file represents a single sample, and peptide-level intensity values are aggregated across all samples to form a consistent dataset suitable for downstream statistical analysis.

To preserve feature-level specificity and enable downstream validation, each peptide entry is assigned a composite identifier (Peptide_ID) derived from protein annotation, m/z value, peptide sequence, modification, and adduct information. This design retains potential ambiguity in protein assignments, allowing later resolution using orthogonal validation methods such as LC-MS/MS.

Input:
- Region-specific Peptide_Summary.csv files (from Step 1)

Process:
1. Read each CSV file as one sample
2. Validate required columns
3. Construct unique Peptide_ID per feature
4. Convert long-format data → wide matrix
5. Fill missing values with 0
6. Sort rows and columns for consistency

Output:
- [Region]_Raw_matrix.csv
- Skipped file report (if any)

❗**Key assumptions (must be understood before use)**
This pipeline assumes that each input file contains the following columns:
```python
{"Protein", "mz", "Peptide", "Modification", "adduct", "Intensity"}
```
This script is not plug-and-play. The following elements must be adjusted depending on your project and this script is written for a single region, other regions requires the adjustment of folder or source files:
```python
data_dir = Path(r"YOUR_PATH_TO_REGION_FOLDER")
output_file = "Region_Raw_matrix.csv" #Region in this project is GM, WM or LC
```

## Step 4 – Data Thresholding and Quality Control
[Data Thresholding and Quality Control.py](https://github.com/user-attachments/files/27372205/Data.Thresholding.and.Quality.Control.py)

This step performs quality control (QC) and filtering on the reconstructed peptide-by-sample matrix to remove low-quality samples and low-confidence features. The goal is to improve data robustness while preserving biologically meaningful signal for downstream statistical analysis.

The QC procedure is applied sequentially: sample-level filtering is performed first to remove poorly detected samples, followed by peptide-level filtering based on detection frequency across the remaining dataset. This order ensures that feature filtering is not biased by low-quality samples.

Input:
- [Region]_Raw_matrix.csv (from Step 2)

Process:
1. Validate matrix structure and numeric values
2. Define detection (presence) based on intensity > 0
3. Filter samples (minimum detected peptides)
4. Filter peptides (minimum detection frequency)
5. Rename sample identifiers (standardized format)
6. Export final QC matrix and reports

Output:
- FINAL_[Region]_QC_matrix.csv
- Removed samples list
- Removed peptides list
- QC summary reports

❗**Key assumptions (must be understood before use)**

1. Detection definition
```python
presence_threshold = 0.0
```
A peptide is considered detected 
- if: **intensity > 0**

This assumption is valid for SCiLS Lab TIC-normalized MALDI-MSI data Zero represents non-detection or below detection limit

2. Sample-level filtering
```python
min_sample_hits = 20
```
Samples are retained only if they contain: **≥ 20 detected peptides** 

3. Peptide-level filtering (core threshold)
   ```python
   min_fraction = 0.40
   min_peptide_samples = ceil(0.40 × number_of_samples)
   ```
A peptide is retained 
- if detected in: **≥ 40% of remaining samples**

This threshold is: empirically determined. Chosen to balance: feature retention and statistical robustness. It is not universal and should be re-evaluated for other datasets.

**4. Sample renaming (project-specific logic)**

- Original format: 1668_finished_1_1_WM_...
- Converted to: 1668_18_WM

This renaming logic is: specific to this project. If using a different dataset: you must modify the regex pattern and/or the coordinate mapping logic. This script is not plug-and-play. The following elements must be adjusted depending on your project and this script is written for a single region, other regions requires the adjustment of folder or source files:
```python
matrix_file = Path(r"YOUR_PATH_TO_RAW_MATRIX")
#one region at a time (e.g., GM, WM, LC)
```
## Step 5 – Metadata Matching and Integration
This step aligns the QC-filtered peptide matrix with sample-level metadata by matching each matrix column to its corresponding biological annotation. The matching is performed using structured identifiers extracted from sample column names, ensuring that each sample is correctly linked to metadata such as tissue diagnosis and experimental variables.

To support region-specific analysis, metadata entries are expanded to include region-level identifiers, allowing direct mapping between matrix columns and anatomical regions (GM, WM, LC). The final output is a metadata table that is fully synchronized with the matrix and ready for downstream statistical modeling. Those identifiers are specific for this project, an only applied if 

Input:
- FINAL_[Region]_QC_matrix.csv (from Step 4)
- Metadata file (FinalDatabase.csv)

Process:
- Validate metadata structure and required columns
- Parse matrix column names into (TMA, Core, Region)
- Expand metadata to region-level entries
- Match matrix samples to metadata
- Identify unmatched or removed samples

Output:
- metadata_MATCHED.csv
- metadata_removed_rows.csv

❗**Key assumptions (must be understood before use)**
1. Matrix column naming format

Matrix columns must follow:[TMA]__[Core]__[Region]. Example: 1668_18_GM.
This format is required for parsing:
```python
r"^(?P<tma>\d+)_(?P<core>\d+)_(?P<region>GM|WM|LC)$"
```

2. Required metadata columns
Metadata must contain: ["TMA", "Core", "TissueDiagnosisID"]

3. TissueRegionID construction
The TissueRegionID construction used in this project is specific to this dataset, where each regional identifier is required to be unique and non-repetitive to ensure accurate mapping between matrix columns and metadata entries.
```python
TissueRegionID = TissueDiagnosisID + "_" + Region
```
This assumes:
- Each core contains all regions (GM, WM, LC)
- Region labels match those in the matrix

This creates a unique identifier for region-level analysis.

4. Working directory
```python
os.chdir("YOUR_SUMMARY_FOLDER_PATH")
matrix_path = "FINAL_GM_QC_matrix.csv"
meta_path = "FinalDatabase.csv"
```
Modify if:
- using WM or LC datasets
- metadata file name differs

## Step 6 (Optional) – Global Feature Extraction for Validation
[Matching Metadata.py](https://github.com/user-attachments/files/27377982/Matching.Metadata.py)

This optional step generates a comprehensive list of unique peptide features across all HiTMaP output files, independent of region-specific processing. Unlike the main matrix construction workflow. This step is designed to support feature-level validation and cross-platform comparison, rather than statistical modeling. By retaining all possible feature combinations and explicitly reporting ambiguous mappings, it enables more rigorous downstream validation against independent datasets such as LC-MS/MS.

The inclusion of mz_align expands feature resolution and ensures compatibility with alignment-based validation workflows.

The pipeline aggregates all input files into a unified dataset and defines feature uniqueness based on a composite combination of protein annotation, m/z value, peptide sequence, modification, adduct, and aligned m/z (mz_align). By incorporating both measured and aligned m/z values, this step increases feature granularity and enables more exhaustive comparison with orthogonal validation data such as LC-MS/MS.

In addition, the pipeline explicitly identifies cases where a single m/z value is associated with multiple peptide sequences, generating a conflict report to highlight ambiguous or potentially misassigned features.

Input:
- All HiTMaP output CSV files (across all regions)

Process:
- Read and validate required columns
- Aggregate all files into a master dataset
- Clean and standardize feature fields
- Identify m/z-to-peptide conflicts
- Define unique features using extended feature definition
- Sort and export full feature list

Output:
- ALL_UNIQUE_FEATURES_LIST.csv
- MZ_PEPTIDE_CONFLICT_REPORT.txt
- SKIPPED_FILES_REPORT.txt

❗**Key assumptions (must be understood before use)**
1. Extended feature definition
```python
#required column:
["Protein", "mz", "Peptide", "mz_align", "desc", "Modification", "adduct"]
#which:
Feature = Protein + mz + Peptide + Modification + adduct + mz_align + desc
```
Compared to the main pipeline:
- includes mz_align
- includes desc
- retains more granular distinctions

This will produce:
- more features
- higher redundancy
- increased ambiguity (intentional for validation)

2. Input directory
```python
in_dir = Path(r"YOUR_PATH_TO_ALL_FILES")
```
