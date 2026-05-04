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
- Recursively search all subfolders for Peptide_Summary.csv files
- Extract core and sample identifiers from folder hierarchy
- Rename files into a standardized format
- Copy all files into a single destination folde

the code loop inside folder that been setted up by HitMaP:
```python
summary_folder = peptide_file.parent
sample_folder = summary_folder.parent
core_folder = sample_folder.parent
```
for example:
```python
E:\Finished Core
├── 1173 finished
│   ├── Sample A
│   │   └── Summary
│   │       └── Peptide_Summary.csv
│   ├── Sample B
│       └── Summary
│           └── Peptide_Summary.csv
├── 1174 finished
    ├── Sample C
        └── Summary
            └── Peptide_Summary.csv
......
```


## Step 2- 
[Separating Files into Brain Regions.py](https://github.com/user-attachments/files/27374390/Separating.Files.into.Brain.Regions.py)


This step organizes raw HiTMaP output files (Peptide_Summary.csv) into region-specific folders based on brain region annotations embedded in file names. All files are initially placed in a single directory and then automatically sorted into Gray Matter (GM), White Matter (WM), and Locus Coeruleus (LC) subfolders.

The script scans each file name for region identifiers and copies the files into corresponding directories, creating a structured dataset for downstream region-specific analysis.

```python
_GM_
_WM_
_LC_
```
Your file has to be named with those region identifiers or you can alter it for your project 

Input:
- Folder containing all HiTMaP Peptide_Summary.csv files

Process:
- Detect region from filename
- Create region-specific folders (GM, WM, LC)
- Copy files into corresponding folders

Output:
ByRegion
- Gray Matter files
- White Matter files
- Locus Coeruleus files

Required user modifications

Before running the script, you must update the file path:
```python
data_dir = Path(r"YOUR_PATH_TO_PEPTIDE_SUMMARY_FILES")

#Example:

data_dir = Path(r"D:\LC_project\HiTMaP_outputs")
```

## Step 2 – Matrix Construction from HiTMaP Output
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

This pipeline assumes that each input file contains the following columns:
```python
{"Protein", "mz", "Peptide", "Modification", "adduct", "Intensity"}
```
This script is not plug-and-play. The following elements must be adjusted depending on your project and this script is written for a single region, other regions requires the adjustment of folder or source files:
```python
data_dir = Path(r"YOUR_PATH_TO_REGION_FOLDER")
output_file = "Region_Raw_matrix.csv" #Region in this project is GM, WM or LC
```

If a comprehensive feature set without region-specific separation is required (e.g., for global feature enumeration or cross-region comparison), please refer to [Step X – Global Feature Aggregation](#step-x--global-feature-aggregation), where all HiTMaP-derived features are consolidated into a unified dataset prior to region-specific processing.

## Step 3 – Data Thresholding and Quality Control
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

**1. Detection definition**
```python
presence_threshold = 0.0
```
A peptide is considered detected 
- if: **intensity > 0**

This assumption is valid for SCiLS Lab TIC-normalized MALDI-MSI data Zero represents non-detection or below detection limit

**2. Sample-level filtering**
```python
min_sample_hits = 20
```
Samples are retained only if they contain: **≥ 20 detected peptides** 

**3. Peptide-level filtering (core threshold)**
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

