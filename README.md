# LComics
Code for data analysis of LC proteomic and glycomic data.

This project investigates region-specific molecular features in the human locus coeruleus (LC) using spatial proteomics and quantitative imaging within a human cohort of 126 individuals without clinically diagnosed neurodegenerative disease. Each tissue microarray (TMA) core was histologically annotated and manually segmented into gray matter (GM), white matter (WM), and LC neuron–enriched regions using QuPath, followed by spatial registration to MALDI mass spectrometry imaging (MALDI-MSI) data processed in SCiLS Lab.

Region-specific mass spectrometry data were exported as normalized spectral files (imzML and ibd) under total ion count normalization, enabling quantitative comparison of molecular signals across spatially defined tissue compartments. Peptide features were subsequently identified using the **HiTMaP** ([GitHub](https://github.com/MASHUOA/HiTMaP) | [Paper](https://doi.org/10.1038/s41467-021-23461-w)) computational framework, which maps MALDI-MSI spectra to in silico peptide libraries while controlling for false discovery.

**This repository implements a downstream computational workflow following HiTMaP output, including large-scale matrix reconstruction, quality filtering, and metadata integration**. The pipeline consolidates region-resolved peptide features into structured matrices, applies empirically defined thresholds to improve robustness, and aligns molecular data with sample-level metadata containing neuropathological and demographic variables such as Braak stage, Thal phase, and APOE genotype.

Together, this workflow provides a reproducible framework for transforming spatial proteomics outputs into analysis-ready datasets, enabling systematic investigation of region-specific molecular alterations associated with early pathological processes in the human brain.
