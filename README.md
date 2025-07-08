# SARS-CoV-2 Genome Diversity Analysis
![Python Version](https://img.shields.io/badge/python-3.7+-blue.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15825646.svg)](https://doi.org/10.5281/zenodo.15825646)

## Project Overview

This repository contains the data, scripts, and results for a large-scale analysis of SARS-CoV-2 nucleotide diversity (π). The primary goal is to provide a reproducible pipeline and a curated dataset for studying viral evolution, particularly the dynamics of selective sweeps and variant transitions.

The core of this project is a validated bioinformatic pipeline that processes the public UShER mutation-annotated phylogeny to generate a whole-genome alignment suitable for population genetic analysis. A comprehensive explanation of the methodology is available in the [Technical Report](./reports/technical_report.md).

## Getting Started

The primary dataset is a 7 GB multiple sequence alignment file that has been split into smaller, compressed parts for ease of access. To reassemble the full dataset, first navigate to the `data/` directory and then use the following command:

```bash
# Ensure you are in the 'data' directory before running
gunzip -c *.fas.gz > all_clades_aligned_combined.fas
```

## Repository Structure

This repository is organized into the following directories:

*   **/data**: Contains the primary compressed dataset files. See the [data/README.md](./data/README.md) for a detailed description of the data, its format, and provenance.
*   **/scripts**: Contains the Python scripts used for the data generation and analysis pipeline. See the [scripts/README.md](./scripts/README.md) for details on dependencies and usage.
*   **/results**: Contains the final output files from the analysis, including key figures and the nucleotide diversity data tables. These are organized into subdirectories based on the analysis (e.g., `jn1_sweep_analysis`, `historical_diversity_analysis`). See the [results/README.md](./results/README.md) for a summary and interpretation of the findings.
*   **/reports**: Contains the detailed `technical_report.md` as well as specific brief reports on various analyses.
    *   `technical_report.md`: A comprehensive document detailing the full data generation pipeline, theoretical background, and analytical methods.
    *   `jn1_sweep_analysis/`: Subdirectory containing brief reports related to the JN.1 lineage analysis:
        *   [Brief Report 1: Quantifying the JN.1 Selective Sweep](./reports/jn1_sweep_analysis/brief_report_1.md)
        *   [Brief Report 2: Comparative Genomics of JN.1 and BA.2.86](./reports/jn1_sweep_analysis/brief_report_2.md)
        *   [Brief Report 3: Visualizing Two Sweeps in JN.1](./reports/jn1_sweep_analysis/brief_report_3.md)
        *   [Brief Report 4: Sampling Bias in F456L Sweep Analysis](./reports/jn1_sweep_analysis/brief_report_4.md)

## Note on Data Processing

For the **historical diversity analysis (2020-2022)**, the nucleotide diversity (π) values in the data file (`results/historical_diversity_analysis/nucdiv_heatmap_data.csv`) were generated with a preliminary version of the analysis script (`scripts/nucdiv_stats_v1_archive.py`). A latent bug in this older script version could cause an artificial inflation of π in the presence of sequences with a high number of ambiguous bases ('N'). The visualization (`results/historical_diversity_analysis/heatmap_nucdiv.png`) presented for this historical period was generated **after** applying a strict quality filter (π < 0.001) to the data, which removed potentially affected data points. All current and future analyses, including the JN.1 study, use a new, fully validated function (`scripts/nucdiv_stats.py` or `scripts/nucdiv_stats_v2.py`) that is robust to this issue.

## Citation

If using this dataset, code, or the methodology in your research, please cite this repository. The recommended citation format is:

*   Friedman, R. (2025). *A Pipeline for Scalable Analysis of SARS-CoV-2 Nucleotide Diversity* (Version 1.0) [Data set and software]. Zenodo. https://doi.org/10.5281/zenodo.15815025

Please also refer to `technical_report.md` for full details on citing the underlying data sources and tools (e.g., UShER, bcftools).

## Acknowledgements

The conceptual development of the methodology and the drafting of the technical and brief reports benefited significantly from discussions and iterative refinement with an AI language model, Gemini 2.5 Pro (Google). The author oversaw and reviewed the accuracy and robustness of all parts of this study.
