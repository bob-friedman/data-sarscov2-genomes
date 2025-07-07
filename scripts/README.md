# Scripts Directory

This directory contains the Python scripts used to generate the data and perform the nucleotide diversity analyses. The scripts are designed to be run in a two-part workflow.

## Scripts

### Part I: Data Generation

*   **`gemini_jn1.py`** (or `nucdiv_v4.py`)
    *   **Purpose:** The main data generation pipeline. It takes the raw UShER protobuf file as input and performs subsampling and consensus sequence generation.
    *   **Output:** The primary `all_clades_aligned.fas` multiple sequence alignment file.

---

### Part II: Nucleotide Diversity Analysis

#### **Current Version (Recommended for All New Work)**

*   **`nucdiv_stats.py`** (or `nucdiv_stats_v2.py`)
    *   **Status:** **Current, Validated Version.**
    *   **Purpose:** This is the robust, fully validated script for analyzing the alignment file. It takes the large alignment as input, converts it to a high-performance binary format, and calculates nucleotide diversity (π) for each lineage, partitioned into time bins.
    *   **Methodology:** This script uses a per-site calculation method that correctly handles missing data and has been validated against a ground-truth, brute-force calculation. It is the recommended script for all new analyses.

#### **Archived Version (For Historical Reproducibility Only)**

*   **`nucdiv_stats_v1_archive.py`**
    *   **Status:** **Archived, Deprecated.**
    *   **Purpose:** This is the preliminary version of the script used for the historical heatmap analysis (2020-2022). It is retained in this repository for the sole purpose of ensuring full reproducibility of that specific result.
    *   **Important Note:** This version contains a known bug that can cause an artificial inflation of π values in the presence of sequences with a high number of ambiguous bases ('N'). **It should not be used for any new work.**

## Dependencies

The core software dependencies for this pipeline are:

| Tool         | Version | Conda Channel |
| :----------- | :------ | :------------ |
| `usher`      | 0.5.8+  | `bioconda`    |
| `bcftools`   | 1.18+   | `bioconda`    |
| `pandas`     | 2.1.4+  | `conda-forge` |
| `biopython`  | 1.81+   | `conda-forge` |
| `scikit-allel`| 1.3.7+  | `conda-forge` |

## Usage

The recommended workflow is a two-step process:
1.  Run the **`gemini_jn1.py`** script to generate the alignment.
2.  Run the current, validated **`nucdiv_stats.py`** script on the output of step 1 to calculate diversity.

For a complete breakdown of the code and methodology, please refer to the main [Technical Report](../technical_report.md).
