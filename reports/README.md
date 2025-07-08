# Comparative Genomic and Temporal Analysis of SARS-CoV-2 Lineage JN.1

## Overview

This repository includes the computational pipeline and findings from a genomic analysis of the SARS-CoV-2 JN.1 lineage. The project investigates the evolutionary divergence of JN.1 from its parental lineage, BA.2.86, to identify genetic adaptations that may explain its rapid global spread.

The above analysis is divided into two main parts, corresponding to two key research questions:

1.  **Comparative Genomics:** How does the genomic diversity of JN.1 differ from BA.2.86? This is addressed through a genome-wide sliding window analysis of nucleotide diversity (π) and a targeted scan for differentiating non-synonymous mutations.

2.  **Temporal Dynamics:** How can we visualize evolutionary events in real-time? This is addressed by tracking the frequency of two key mutations over time within the JN.1 lineage: the lineage-defining `L455S` mutation and the emerging `F456L` mutation.

The findings from these analyses are detailed in the `technical_report.md` located in this directory, and further summarized in specific brief reports.

## Directory Structure

This `reports/` directory is organized as follows:

*   **`technical_report.md`**: A comprehensive document detailing the full data generation pipeline, theoretical background, and analytical methods used across the project.
*   **`jn1_sweep_analysis/`**: This subdirectory contains a collection of brief reports focusing on various aspects of the SARS-CoV-2 JN.1 lineage analysis, including its selective sweep, comparative genomics with BA.2.86, and investigations into specific mutations.
    *   `brief_report_1.md`: Quantifying the Selective Sweep of the SARS-CoV-2 JN.1 Lineage.
    *   `brief_report_2.md`: A Comparative Genomic Analysis of SARS-CoV-2 Lineages JN.1 and BA.2.86.
    *   `brief_report_3.md`: Visualizing Viral Evolution in Real-Time: A Tale of Two Sweeps in the SARS-CoV-2 JN.1 Lineage.
    *   `brief_report_4.md`: Distinguishing a Localized Outbreak from a Global Selective Sweep: A Case Study on Sampling Bias.
*   **`README.md`**: This file, providing an overview of the reports directory.

For details on the scripts used to generate these analyses and the results/plots, please refer to the `scripts/` and `results/` directories in the main repository, respectively.
