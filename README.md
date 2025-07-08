[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15825646.svg)](https://doi.org/10.5281/zenodo.15825646)
# SARS-CoV-2 Genome Diversity Analysis

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
*   **/results**: Contains the final output files from the analysis, including key figures and the nucleotide diversity data tables. See the [results/README.md](./results/README.md) for a summary and interpretation of the findings.
*   **technical_report.md**: A comprehensive document detailing the full data generation pipeline, theoretical background, and analytical methods.
*   **brief_report_1.md**: A [Brief Report](./reports/brief_report_1.md) that describes the quantification of the selective sweep event of the JN.1 lineage in the SARS-CoV-2 population.
*   **brief_report_2.md**: A [Brief Report](./reports/brief_report_2.md) that describes a comparative genomic analysis in the context of adapatation in the lineages JN.1 and BA.2.86 in the SARS-CoV-2 population.
*   **brief_report_3.md**: A [Brief Report](./reports/brief_report_3.md) that relies on a codon analysis to characterize two distinct selective sweeps relevant to the SARS-CoV-2 JN.1 lineage.
*   **brief_report_4.md**: A [Brief Report](./reports/brief_report_4.md) shows how sampling heterogeneity can lead to analytical artifacts in the characterization of a selective sweep in the SARS-CoV-2 JN.1 lineage.

## Note on Data Processing

The nucleotide diversity (π) values in the initial raw results file (`results/nucleotide_diversity_results_filtered.csv`) were generated with a preliminary version of the analysis script. A latent bug in this version could cause an artificial inflation of π in the presence of sequences with a high number of ambiguous bases ('N'). The primary visualization (heatmap) was generated **after** applying a strict quality filter (π < 0.001), which removed all affected data points. All current and future analyses use a new, fully validated function that is robust to this issue.

## Citation

If using this dataset, code, or the methodology in your research, please cite this repository. The recommended citation format is:

*   Friedman, R. (2025). *A Pipeline for Scalable Analysis of SARS-CoV-2 Nucleotide Diversity* (Version 1.0) [Data set and software]. Zenodo. https://doi.org/10.5281/zenodo.15815025

Please also refer to `technical_report.md` for full details on citing the underlying data sources and tools (e.g., UShER, bcftools).

## Acknowledgements

The conceptual development of the methodology and the drafting of the technical and brief reports benefited significantly from discussions and iterative refinement with an AI language model, Gemini 2.5 Pro (Google). The author oversaw and reviewed the accuracy and robustness of all parts of this study.
