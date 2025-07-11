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

---

### Standalone Analysis Script (for specific JN.1 reports)

*   **`analysis_pipeline.py`**
    *   **Purpose:** This script performs a specific series of analyses related to the initial JN.1 vs BA.2.86 comparison and temporal tracking of L455S/F456L mutations, as detailed in `brief_report_2.md` and `brief_report_3.md`. It handles data loading (expecting `all_clades_aligned.fas` and other files to be present), diversity calculations, and plot generation for these specific reports.
    *   **Note:** While this script is functional for its specific purpose, the combination of `gemini_jn1.py` and `nucdiv_stats.py` is recommended for newer, more flexible analyses.

## Usage

The general recommended workflow is a two-step process:
1.  Run the **`gemini_jn1.py`** script (or `nucdiv_v4.py`) to generate the `all_clades_aligned.fas` alignment.
2.  Run the current, validated **`nucdiv_stats.py`** script (or `nucdiv_stats_v2.py`) on the output of step 1 to calculate nucleotide diversity for various lineages and time periods.

Alternatively, to reproduce the specific analyses for the JN.1 vs BA.2.86 comparison and L455S/F456L temporal trends (as in `brief_report_2.md` and `brief_report_3.md`), you can run **`analysis_pipeline.py`**. Ensure all prerequisite files listed in its header are available in the same directory.

For a complete breakdown of the code and methodology, please refer to the main [Technical Report](../reports/technical_report.md).

---

### Standalone Script: Ancestral Sequence Reconstruction

*   **`parsimony_reconstruction.py`**
    *   **Purpose:** This script reconstructs ancestral biological sequences on a phylogenetic tree using a parsimony-based method. It takes a Newick tree file with mutational changes annotated on branches and a reference sequence for the root node as input.
    *   **Input:** A text file where the first line is the Newick string (e.g., `((Leaf1:[1C>A],Leaf2:[1C>T])Anc1)Root;`) and the second line is the root/reference sequence (e.g., `CGTGA`). Mutations are 0-indexed.
    *   **Output:** Prints the inferred sequences for all internal nodes and leaves to standard output.
    *   **Methodology:**
        1.  **Initial Pass (Root to Tips):** Calculates sequences for all nodes by applying specified mutations down each branch, starting with the provided root sequence.
        2.  **Second Pass (Tips to Root):** Refines internal node sequences using Fitch's parsimony algorithm. For each site, if the sets of possible nucleotides from children nodes intersect, the intersection is assigned. If not, the union is assigned. IUPAC ambiguity codes are used to represent sets of nucleotides.
    *   **Usage:** `python scripts/parsimony_reconstruction.py <path_to_input_file>`
