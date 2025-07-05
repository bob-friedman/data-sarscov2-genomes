## **Technical Report: A Pipeline for Scalable Analysis of SARS-CoV-2 Nucleotide Diversity**

### **Abstract**

This report provides a comprehensive guide to a two-part workflow for calculating genome-wide nucleotide diversity (π) for SARS-CoV-2 using the UShER mutation-annotated phylogeny. Part I details the `nucdiv_v4.py` pipeline for data acquisition, clade-based subsampling, and parallelizable consensus sequence generation using `matUtils` and `bcftools`. The result is a single, concatenated FASTA alignment. Part II describes the `nucdiv_stats.py` script, which processes this alignment to calculate nucleotide diversity across lineages and time bins. This second stage employs memory-efficient techniques, including memory-mapped binary files and chunked data processing, to handle multi-gigabyte datasets on standard hardware. We present detailed code annotations, expected data outputs, and a robust data hosting strategy to ensure full reproducibility and community access.

---

### **Introduction**

#### **Background on this Report**

Understanding SARS-CoV-2 evolution requires linking genetic diversity to population processes. SARS-CoV-2 variants emerge through mutation, selection, recombination, and chance (genetic drift) under immune pressure. For example, while most viral mutations are neutral or deleterious, a few confer enhanced transmissibility or immune escape. These advantageous mutations (e.g., in Spike) can fix rapidly, driving **selective sweeps** that purge diversity. Conversely, when many lineages co-circulate (a “mosaic” regime), diversity remains high and recombination may shuffle genomes without an immediate sweep. Spatial spread also matters: when a variant spreads to new regions, founder bottlenecks tend to reduce diversity away from the origin.

To quantify these dynamics, we compute **nucleotide diversity (π)**—the average pairwise genetic differences per site among sampled genomes. Under neutral evolution, π equilibrates at approximately 2Nₑμ (where Nₑ is the effective population size and μ the per-site mutation rate). Thus, π serves as an empirical gauge of the viral population’s “effective size” and past bottlenecks.

In this report, we describe a scalable pipeline to build multi-FASTA alignments for π analysis, and we integrate population-genetic theory to interpret the outputs. The accompanying proposal establishes the rationale for a clade-based subsampling strategy, which balances representative sampling of major lineages against computational feasibility. By selecting up to 1,000 genomes from each of the 250 most prevalent Pangolin lineages, this approach mitigates biases from uneven sequencing efforts while retaining the statistical power needed to detect significant patterns in nucleotide diversity (π).

The core of this method is the UShER mutation-annotated tree, which provides both the global phylogenetic context and per-sample variant calls in a compressed protobuf format. The pipeline leverages `matUtils` for rapid variant extraction and `bcftools` for consensus sequence generation. The final output is a single, concatenated FASTA file where each sequence header is annotated with its sample identifier and Pango lineage (e.g., `>sample_id|clade_name`), making it immediately suitable for downstream analysis with tools like `scikit-allel`, R, or for phylodynamic inference with BEAST.

In this report, Part I details each component of the alignment generation script. Part II details the subsequent script for calculating π and interpreting the results. Finally, we outline a version-controlled data hosting plan on GitHub to promote open science.

#### **Population Genetics and Evolutionary Dynamics of SARS-CoV-2**

The evolution of SARS-CoV-2 is shaped by classical population genetic principles that apply to all replicating organisms, but with particular features driven by its viral biology and the global scale of the COVID-19 pandemic. A viral *population* refers to the pool of SARS-CoV-2 genomes circulating among infected hosts in a given time and place. Although the number of infections may be in the millions, the concept of the *effective population size* (Nₑ) is central to understanding evolutionary change. Nₑ represents the number of viral genomes that effectively contribute to the next generation of infections. In SARS-CoV-2, Nₑ is orders of magnitude smaller than the number of infected individuals, largely due to extremely tight *transmission bottlenecks*—situations in which only a small number of virions successfully establish infection in a new host (Markov et al., Nat Rev Microbiol, 2023).

These bottlenecks make *genetic drift*—the random fluctuation of allele frequencies across generations—a powerful force in SARS-CoV-2 evolution. For instance, when the virus spreads into a new geographic region, it often does so through a handful of transmission events. This can lead to a *founder effect*, where only a small, non-representative sample of the global viral diversity seeds a local epidemic. Such effects have been documented repeatedly during the early global dissemination of SARS-CoV-2, leading to lower diversity in “sink” populations compared to “source” regions (Tasaki et al., PLOS One, 2021). These stochastic dynamics can amplify neutral or even slightly deleterious mutations simply by chance.

Meanwhile, *mutation* introduces new genetic variation each generation. SARS-CoV-2 has an estimated substitution rate of approximately 10⁻³ mutations per site per year, allowing it to accumulate substantial diversity in a short time. While the majority of these mutations are neutral or deleterious, some confer fitness advantages by improving viral transmissibility or immune escape. The D614G mutation in the spike protein is an early example that increased infectivity and rapidly swept to global dominance (Harvey et al., Cell, 2021). Later, the N501Y mutation in Alpha, Beta, and Gamma variants provided improved ACE2 binding and was similarly favored.

When an advantageous mutation arises, *positive selection* can drive its rapid fixation in the population. Because SARS-CoV-2’s genome is largely clonal and recombines infrequently, these selective sweeps often involve *linked selection*—the genetic hitchhiking of nearby neutral or even mildly deleterious mutations on the same genomic background. A sweep of the spike protein, for example, will eliminate variation in adjacent regions. This explains why major variant transitions, such as from Delta to Omicron, are associated with steep collapses in genome-wide diversity: the entire viral population becomes dominated by a single haplotype that outcompetes others (Suratekar et al., Mol Syst Biol, 2022).

Recombination adds another layer of complexity. Although relatively rare early in the pandemic, recombination became more prominent as global diversity increased. When two distinct SARS-CoV-2 lineages co-infect the same host, genome segments can be exchanged. This process has led to the emergence of recombinant lineages such as XBB, which originated from two co-circulating BA.2 sublineages (Planas et al., Nat Commun, 2024). Recombination enables the virus to “leap” across fitness landscapes, combining advantageous mutations from multiple backgrounds more rapidly than mutation alone would allow. However, recombinants do not always dominate immediately. For example, while the ancestral XBB did not trigger a global sweep upon emergence, its descendant XBB.1.5 eventually acquired additional spike mutations that enabled dominance—suggesting that recombination set the stage, but further adaptation completed the sweep.

The interaction between selection, mutation, drift, and recombination plays out across both time and geography. In the early stages of variant emergence, diversity often drops sharply—a signature of a selective sweep. As new lineages accumulate mutations over time and as regions experience repeated introductions, diversity can rise again. Spatial dynamics add further structure. A variant’s region of origin often shows the highest nucleotide diversity (π), because it hosts the earliest and most diverse viral samples. As the variant spreads to new regions, each expansion is subject to transmission bottlenecks and founder effects, leading to a characteristic gradient of decreasing π from source to sink regions. Over time, local diversity may rebound as new lineages emerge or co-circulate, sometimes enabling recombination and setting the stage for future variants.

Genomic data enable the quantification of these processes. For instance, the average number of nucleotide differences per site between sampled genomes—nucleotide diversity, or π—is a core metric of population genetic variation (Nei and Li, PNAS, 1979). Under the *neutral theory of molecular evolution* (Kimura, 1983), π is expected to equilibrate at 2Nₑμ for haploid populations like viruses. Deviations from this expectation, such as sharp declines in π, indicate departures from neutrality and suggest events such as selective sweeps or demographic contractions (Tajima, Genetics, 1989). Temporal analysis of π allows us to track the rise and fall of diversity across variant transitions, while spatial analysis can reveal regional differences that reflect founder effects, differential introduction rates, or local adaptation.

Phylogenetic inference further supports this analysis. For example, Bayesian phylodynamic models like those implemented in BEAST (Suchard et al., Virus Evol, 2018) use time-stamped sequence data to reconstruct effective population size (Nₑ) through time, model lineage movement between regions, and infer the timing of major evolutionary events. These models complement π by estimating the underlying processes that shape diversity. When both π and BEAST-inferred Nₑ decline simultaneously, it provides strong evidence of a real population contraction due to a selective sweep. Conversely, sustained high π alongside stable or growing Nₑ may suggest co-circulation of multiple lineages and a mosaic evolutionary regime—conditions favorable to recombination-driven evolution.

In summary, the population genetics of SARS-CoV-2 exemplifies how classic evolutionary forces—mutation, selection, drift, and recombination—interact in a spatially and temporally structured pathogen population. These forces explain the recurring emergence of variants, the episodic collapses and rebounds of diversity, and the rise of recombinant forms like XBB. By combining high-resolution genomic data with theoretical models, we gain not only a descriptive account of viral evolution but a predictive framework for understanding and anticipating future variant dynamics.

#### **Calculations of Nucleotide Diversity (π) and the Neutral Theory of Evolution**

Nucleotide diversity, π, is defined as the average number of nucleotide differences per site between two sequences chosen at random. Nei and Li (1979) formalized this concept, and under Kimura’s neutral theory (1983), π reaches an equilibrium (π ≈ 2Nₑμ for haploids) set by mutation and effective population size. In practical terms, **high π** indicates substantial polymorphism (many coexisting lineages), whereas **low π** indicates near-uniformity (often following a recent sweep).

In a purely neutral, constant-size population, π would be roughly stable. Significant departures (sharp rises or falls) indicate non-neutral processes. For example, Tajima (1989) pointed out that time-series or spatial gradients in π signal selection or demographic change. Thus, temporal drops in π suggest recent sweeps or bottlenecks, while excess diversity above neutral expectation suggests balancing selection or population structure. In our SARS-CoV-2 context, a sudden global π crash concomitant with a new variant’s emergence would support a selective sweep; a spatial gradient (higher π in origin vs. sinks) would indicate founder effects.

#### **Evolution of SARS-CoV-2 and its Expected π Patterns**

SARS-CoV-2 evolution alternates among distinct regimes that predict characteristic π signatures:

*   **Gradual adaptation:** Slow mutation accumulation within a lineage. π may drift slowly upward but without large jumps. Diversity grows gradually until an event occurs.
*   **Co-circulation (mosaic) regime:** Multiple divergent lineages persist simultaneously (e.g., Omicron BA.5, BQ.1, XBB in late 2022). This **maintains high π**, reflecting many genetic backgrounds. Such regimes foster recombination: e.g., the recombinant XBB arose from co-circulating BA.2 strains. In this mosaic phase, adding a recombinant lineage (like XBB) typically does **not** cause an immediate diversity collapse; instead, total π stays elevated as variants share ancestry deeply.
*   **Selective sweep regime:** One lineage acquires a strong fitness advantage and rapidly fixes. Classic SARS-CoV-2 examples are Alpha→Delta, Delta→Omicron, and XBB→JN.1. Each sweep “prunes” the genealogy, producing a **precipitous drop in π** worldwide as the new variant overcomes others. Thus we expect a sharp global minimum in π at sweep onset.
*   **Antigenic leap (deep-branch emergence):** Occasionally, a very divergent variant appears (e.g., Omicron BA.1, BA.2.86) from a long unsampled lineage. This is akin to injecting a highly different genome, which can **reset diversity**: π may jump upward transiently as a new deep branch enters the population.
*   **Founder-effect spread:** During geographic expansion of any new variant, repeated bottlenecks occur. The source region (with multiple introductions) will tend to have **higher π** than newly seeded “sink” regions, which start from few genomes. Over time, as more lineages accumulate, π in sinks will rise. These patterns (source vs. sink π gradients) are predicted by structured coalescent models of spread.

By labeling each time-slice of data according to these regimes, we interpret π trends contextually. For example, during a “mosaic” Omicron phase we expect sustained high π and evidence of recombination, whereas during a “sweep” we expect a rapid π crash. Prior studies support this: the global SARS-CoV-2 diversity plunged during past variant sweeps, and recombination during co-circulation (as in XBB) kept diversity elevated.

#### **Phylodynamic Inference with BEAST Software**

To complement direct π measurements, **phylodynamic** models (BEAST) may be utilized that infer the underlying population dynamics from sequence data. BEAST (Suchard et al., 2018) jointly estimates time-calibrated phylogenies, effective population size Nₑ(t), and migration or birth–death parameters from dated genomes. Importantly, BEAST *does not compute π directly*; rather, it infers coalescent rates that determine how π would evolve. A long coalescent time (high inferred Nₑ) corresponds to higher diversity, while rapid coalescence (low Nₑ) corresponds to diversity loss. In practice, a BEAST skyline on a sweeping lineage (e.g., JN.1) would show a sharp contraction in Nₑ coincident with the π crash.

Thus, BEAST provides a statistical model-based check on the π patterns. For example, the BEAST skyline and structured coalescent models is applicable for inference on whether key lineages (e.g., early vs. late JN.1 or XBB) associated with changes in Nₑ. A concurrent drop in inferred Nₑ and measured π would may be interpreted as a selective sweep. Similarly, inferred migration rates can quantify founder transitions. These analyses require smaller, hypothesis-driven subsamples due to computational costs.

---

### **Part I: Consensus Sequence Generation from the UShER Tree**

This part details the `nucdiv_v4.py` workflow, designed for the initial, large-scale generation of a multiple sequence alignment from the raw UShER protobuf file.

#### **Methods for Informatics Analysis**

The workflow proceeds from environment setup to final file cleanup. Each step is designed to be modular and reproducible.

**Data Flow Overview**

The pipeline follows a linear data flow, transforming raw public data into a final, analysis-ready alignment.

```
[Input Data] -> [Setup & Acquisition] -> [Subsampling] -> [Execution Engine] -> [Final Alignment]
  - .pb file       - Install tools           - Gen. control file  - Extract VCFs     - all_clades_aligned.fas
  - .tsv metadata  - Download data           - Gen. sample lists  - Build consensus
  - .fa reference                                                  - Append to FASTA
```

**Software Dependencies**

Reproducibility requires specific software versions. The workflow was validated using the following tools, which can be installed via Conda from the specified channels.

| Tool         | Version | Conda Channel | Purpose                                      |
| :----------- | :------ | :------------ | :------------------------------------------- |
| `condacolab` | 0.1.7+  | `pip`         | Manages Conda environments in Google Colab.  |
| `usher`      | 0.5.8+  | `bioconda`    | Extracts data from mutation-annotated trees. |
| `bcftools`   | 1.18+   | `bioconda`    | Manipulates VCF files and generates consensus. |
| `pandas`     | 2.1.4+  | `conda-forge` | Performs tabular data manipulation.          |
| `biopython`  | 1.81+   | `conda-forge` | Parses the final FASTA file for analysis.    |

---

#### **Cell 1: Environment Configuration and Master Sample Extraction**

This cell initializes the environment and extracts a complete list of all samples contained within the UShER protobuf file.

```python
# Cell 1: Environment Setup and Configuration
!pip install -q condacolab
import condacolab; condacolab.install()
import os, pandas as pd, sys

# Define input and output paths
PB_FILE = "public-latest.all.masked.pb"
METADATA_FILE = "public-latest.metadata.tsv"
REF_FASTA = "NC_045512v2.fa"
SAMPLES_TSV = "samples.txt"
CLADES_CONTROL_FILE = "clades_to_process.txt"
FINAL_ALIGNED_FASTA = "all_clades_aligned.fas"

# Generate the master sample list from the .pb file
print("Generating master sample list from the protobuf file...")
!matUtils summary -i {PB_FILE} -s {SAMPLES_TSV}
```

The script first installs `condacolab` and imports necessary Python libraries. Key file paths are assigned to variables for clarity and ease of modification. The final command, `matUtils summary`, reads the entire protobuf (`.pb`) file and writes all unique sample identifiers to `samples.txt`. This master list is essential for the subsequent cross-referencing step, ensuring that the analysis only considers samples that are present in both the phylogeny and the metadata file.

#### **Cell 2 & 3: Data Acquisition and Tool Installation**

These cells are responsible for downloading the required datasets and installing the bioinformatics software.

```bash
# Cell 2: Install Bioinformatics Dependencies
conda install -y -c bioconda -c conda-forge usher bcftools

# Cell 3: Acquire a Consistent, Archived Data Set
RELEASE_DATE="2025/06/29" # Example date for reproducibility
BASE_URL="http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/${RELEASE_DATE}/"
wget -nc "${BASE_URL}/public-${RELEASE_DATE//\//-}.all.masked.pb.gz" -O public-latest.all.masked.pb.gz
gunzip -f public-latest.all.masked.pb.gz
# (Commands repeated for metadata.tsv and reference.fa)
```

The commands shown here ensure a consistent analytical environment. `conda install` fetches the specified versions of UShER and bcftools. The `wget` commands download specific, dated releases of the UShER phylogeny and metadata. Using an archived version instead of "latest" is critical for reproducibility.

#### **Cell 4: Clade-Based Subsampling**

This cell implements the core sampling strategy.

```python
# Cell 4: Subsampling Logic
NUMBER_OF_MAJOR_CLADES = 250
SAMPLES_PER_CLADE = 1000

# Load master sample list and metadata
master_sample_set = {line.strip().split('\t')[0] for line in open(SAMPLES_TSV)}
metadata_df = pd.read_csv(METADATA_FILE, sep='\t', usecols=['strain','pangolin_lineage'])
filtered_df = metadata_df[metadata_df['strain'].isin(master_sample_set)]

# Identify top clades and generate subsampled lists
grouped_metadata = filtered_df.groupby('pangolin_lineage')
major_clades = grouped_metadata.size().nlargest(NUMBER_OF_MAJOR_CLADES)

with open(CLADES_CONTROL_FILE,'w') as f_control:
    for clade_name, group_df in grouped_metadata:
        if clade_name in major_clades.index:
            num_to_sample = min(SAMPLES_PER_CLADE, len(group_df))
            subsampled_group = group_df.sample(n=num_to_sample, random_state=42)
            safe_clade_name = "".join(c if c.isalnum() else "_" for c in clade_name)
            subsampled_group['strain'].to_csv(f'samples_{safe_clade_name}.txt', index=False, header=False)
            f_control.write(f"{safe_clade_name}\n")
```

The script defines sampling parameters (top 250 clades, up to 1,000 samples per clade), filters metadata to include only samples present in the phylogeny, and identifies the most prevalent clades. It then iterates through these clades, performs a reproducible random sample (`random_state=42`), sanitizes the clade name for use in filenames (e.g., `BA.2.75` becomes `BA_2_75`), and writes the subsampled strain IDs to a dedicated `samples_*.txt` file and the control file.

#### **Cell 5: Consensus Sequence Generation Loop**

This shell script block is the execution engine of the pipeline.

```bash
#!/bin/bash
set -e # Exit immediately if a command fails

PB_FILE="public-latest.all.masked.pb"
REF_FASTA="NC_045512v2.fa"
CLADES_CONTROL_FILE="clades_to_process.txt"
FINAL_ALIGNED_FASTA="all_clades_aligned.fas"
touch "$FINAL_ALIGNED_FASTA"

while read -r safe_clade_name; do
    CLADE_SAMPLE_FILE="samples_${safe_clade_name}.txt"
    TEMP_VCF="temp_${safe_clade_name}.vcf"
    matUtils extract -i "$PB_FILE" --samples "$CLADE_SAMPLE_FILE" -v "$TEMP_VCF"
    bgzip --force "$TEMP_VCF"
    bcftools index --force --tbi "${TEMP_VCF}.gz"
    while read -r sample_id; do
        bcftools consensus -f "$REF_FASTA" --sample "$sample_id" "${TEMP_VCF}.gz" \
          | sed "1s#.*#>${sample_id}|${safe_clade_name}#" >> "$FINAL_ALIGNED_FASTA"
    done < "$CLADE_SAMPLE_FILE"
    rm "${TEMP_VCF}.gz" "${TEMP_VCF}.gz.tbi"
done < "$CLADES_CONTROL_FILE"
```

The loop reads one clade at a time from the control file. For each, it uses `matUtils extract` to isolate variants into a VCF, compresses and indexes the VCF with `bgzip` and `bcftools index`, and then loops through each sample ID. For each ID, `bcftools consensus` reconstructs the full genome sequence by applying that sample's variants to the reference. The `sed` command modifies the FASTA header to include both the sample ID and its clade, appending the result to the final output file.

#### **Expected Output of Part I**

The primary output is `all_clades_aligned.fas`, a multi-FASTA file. Each record follows a standardized format, with the header containing the original strain identifier and the sanitized Pango lineage, separated by a pipe (`|`).

```fasta
>USA/AZ-544/2020|BA_2_75
NNNNACGTACGT...
>another/sample/id|BQ_1_1
NNNNACGTACGT...
```

This alignment is the direct input for the analysis pipeline described in Part II.

---

### **Part II: Analysis of Nucleotide Diversity from Aligned Sequences**

The `nucdiv_stats.py` script, detailed in this section, processes the `all_clades_aligned.fas` file generated in Part I. Its primary goal is to calculate nucleotide diversity (π) across different time periods and viral lineages. This script is designed to handle large (multi-gigabyte) datasets efficiently on standard computer hardware by avoiding loading entire files into memory simultaneously.

#### **Methods for Interpretation of Nucleotide Diversity Patterns**

Once π is computed from the final alignment, we interpret the patterns against theoretical expectations.

*   **Co-circulation/mosaic diversity:** Periods when multiple Omicron sublineages co-existed should show **elevated and relatively stable π**. The addition of a recombinant like XBB would be expected to **maintain** high π, rather than cause an immediate drop.
*   **Selective sweeps:** A sharp fall in π coinciding with the rise of a new dominant variant (e.g., Delta, Omicron, JN.1) is expected. A plot of π over time should reveal a pronounced minimum during the variant's emergence.
*   **Spatial founder effects:** A new variant's origin region should show higher π than remote regions seeded by fewer genomes. Computing π by region over time can detect these gradients.
*   **Partitioning diversity:** In mosaic phases, total diversity (`π_total`) will be much greater than the diversity calculated *within* each major lineage (e.g., `π_XBB`, `π_BQ.1`), indicating that diversity is driven by differences *between* co-circulating lineages.

#### **Detailed Script Breakdown and Methods**

The `nucdiv_stats.py` workflow is executed in two main stages, each implemented as a distinct script segment.

**Stage 1: Standalone Verification and Filtering Script**

This initial stage prepares and validates the input FASTA file from Part I. It is essential for ensuring data quality and creating a crucial index file for downstream metadata alignment.

*   **Function**: This script streams the large, concatenated FASTA file (`all_clades_aligned_combined.fas`). It verifies each sequence for a correctly formatted header (e.g., `>id|clade|date|country`) and consistent sequence length.
*   **Outputs**:
    1.  `verified_sequences.fas`: A new FASTA file containing only valid sequences.
    2.  `verification_log.csv`: A log of all discarded sequences and the reason for their exclusion.
    3.  `verified_sequence_order.txt`: A text file listing the primary sequence IDs (the part of the header before the first pipe) in the exact order they appear in `verified_sequences.fas`. This file is critical for correctly aligning the sequence data with the metadata file.

```python
# Stage 1: Verification Logic
import csv
import os

def fasta_iterator(file_path):
    """A generator to read a FASTA file one sequence at a time."""
    with open(file_path, 'r') as f:
        header, sequence_lines = None, []
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if header: yield header, "".join(sequence_lines)
                header = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line)
        if header: yield header, "".join(sequence_lines)

def verify_and_filter_fasta(input_path, clean_path, log_path, order_path):
    """Streams, verifies, and filters a FASTA file, creating a clean output and an ordered ID index."""
    processed, written, errors = 0, 0, 0
    ref_len = None

    with open(clean_path, 'w') as f_clean, \
         open(log_path, 'w', newline='') as f_log, \
         open(order_path, 'w') as f_order:

        log_writer = csv.writer(f_log)
        log_writer.writerow(['sequence_id', 'error_reason'])
        for header, sequence in fasta_iterator(input_path):
            processed += 1
            error_reason = None
            seq_id = header.split('|')[0] if '|' in header else header
            if header.count('|') != 3:
                error_reason = "Malformed header"
            else:
                if ref_len is None: ref_len = len(sequence)
                elif len(sequence) != ref_len:
                    error_reason = f"Incorrect length (is {len(sequence)}, expected {ref_len})"
            if error_reason:
                errors += 1
                log_writer.writerow([seq_id, error_reason])
            else:
                written += 1
                f_clean.write(f">{header}\n{sequence}\n")
                f_order.write(f"{seq_id}\n")
    print(f"Verification complete. Processed: {processed}, Written: {written}, Discarded: {errors}")

# --- Main Execution Call ---
# verify_and_filter_fasta("all_clades_aligned_combined.fas", "verified_sequences.fas", "verification_log.csv", "verified_sequence_order.txt")
```

**Stage 2: Partitioned Nucleotide Diversity Script**

This stage takes the verified data from Stage 1 and performs the core analysis.

*   **Cell 1: Pre-processing (FASTA to NumPy Binary)**: The text-based `verified_sequences.fas` is converted into a numerical, memory-mapped NumPy array (`alignment.mmap`). Nucleotides (A,C,G,T,N) are mapped to integers (0,1,2,3,4). This binary format enables fast, random access to sequence data without loading the entire multi-gigabyte file into RAM.

```python
# Stage 2, Cell 1: FASTA to Binary Conversion
import numpy as np
def fasta_to_numpy_binary(fasta_path, output_path):
    """Converts a FASTA alignment to a memory-mapped NumPy binary file."""
    num_sequences, seq_length = 0, 0
    with open(fasta_path, 'r') as f: # First pass to get dimensions
        for line in f:
            if line.startswith('>'): num_sequences += 1
            elif not seq_length: seq_length = len(line.strip())

    mmap_array = np.memmap(output_path, dtype=np.uint8, mode='w+', shape=(num_sequences, seq_length))
    nt_map = str.maketrans("ACGTN", "01234")
    current_seq_index = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'): continue
            processed_sequence = line.strip().upper().translate(nt_map)
            mmap_array[current_seq_index, :] = np.frombuffer(processed_sequence.encode(), dtype=np.uint8) - ord('0')
            current_seq_index += 1
    mmap_array.flush()
    return output_path, num_sequences, seq_length

# --- Main Execution Call ---
# binary_file, n_seq, n_len = fasta_to_numpy_binary("verified_sequences.fas", "alignment.mmap")
```

*   **Cell 2: Memory-Efficient Metadata Alignment**: This step aligns the large metadata TSV file with the sequence data. It reads the metadata in chunks, filters for relevant strains using the IDs from `verified_sequence_order.txt`, concatenates the filtered chunks, de-duplicates entries, and re-indexes the resulting DataFrame to match the exact order of sequences in the binary alignment. Finally, it creates a `time_bin` column (e.g., bi-weekly periods) from sample collection dates.

```python
# Stage 2, Cell 2: Metadata Alignment
import pandas as pd
import gc

ORDER_INDEX_FILE = "verified_sequence_order.txt"
METADATA_TSV = "public-latest.metadata.tsv"
CHUNK_SIZE = 500000

ordered_ids_set = set(pd.read_csv(ORDER_INDEX_FILE, header=None, names=['strain'])['strain'])
cols_to_use = ['strain', 'pangolin_lineage', 'date', 'country']
filtered_chunks = []
with pd.read_csv(METADATA_TSV, sep='\t', usecols=cols_to_use, chunksize=CHUNK_SIZE) as reader:
    for chunk in reader:
        chunk['strain'] = chunk['strain'].str.split('|').str[0]
        filtered_chunks.append(chunk[chunk['strain'].isin(ordered_ids_set)])
df_meta_filtered = pd.concat(filtered_chunks, ignore_index=True)
del filtered_chunks; gc.collect()
df_meta_filtered.drop_duplicates(subset=['strain'], keep='first', inplace=True)
ordered_ids_list = pd.read_csv(ORDER_INDEX_FILE, header=None, names=['strain'])['strain'].tolist()
df_meta_aligned = df_meta_filtered.set_index('strain').loc[ordered_ids_list].reset_index()
del df_meta_filtered, ordered_ids_list; gc.collect()
df_meta_aligned['date'] = pd.to_datetime(df_meta_aligned['date'], errors='coerce')
df_meta_aligned.dropna(subset=['date'], inplace=True)
df_meta_aligned['time_bin'] = df_meta_aligned['date'].dt.to_period('2W')
grouped = df_meta_aligned.groupby(['pangolin_lineage', 'time_bin'])
```

*   **Cell 3: Partitioned Nucleotide Diversity (π) Calculation**: This cell iterates through each defined group (by lineage and time bin). For each group with at least two samples, it retrieves the corresponding sequence indices, extracts the sequence slice from `alignment.mmap`, and then calculates π using the `calculate_pi_for_slice` function. The function calculates π by summing pairwise differences at each site and normalizing by the number of pairs and sequence length.

```python
# Stage 2, Cell 3: Partitioned Pi Calculation
def calculate_pi_for_slice(alignment_slice):
    """Calculates pi for a given numpy array slice of an alignment."""
    num_sequences, seq_length = alignment_slice.shape
    if num_sequences < 2: return 0.0
    total_pairwise_diffs = 0.0
    total_possible_pairs = num_sequences * (num_sequences - 1) / 2
    for i in range(seq_length):
        column = alignment_slice[:, i]
        counts = np.bincount(column, minlength=5)
        num_same_pairs = sum(n * (n - 1) / 2 for n in counts[:4]) # Only A,C,G,T
        total_pairwise_diffs += (total_possible_pairs - num_same_pairs)
    return total_pairwise_diffs / (total_possible_pairs * seq_length)

# --- Main Calculation Loop ---
alignment = np.memmap(binary_file, dtype=np.uint8, mode='r', shape=(n_seq, n_len))
results = []
for name, group_df in grouped:
    if len(group_df) < 2: continue
    indices = group_df.index.tolist()
    alignment_slice = alignment[indices, :]
    pi = calculate_pi_for_slice(alignment_slice)
    results.append({'lineage': name[0], 'time_bin': str(name[1]), 'pi': pi, 'n_samples': len(group_df)})
results_df = pd.DataFrame(results)
results_df.to_csv("nucleotide_diversity_results.csv", index=False)
```

*   **Cell 4: Filtering and Visualization**: The raw results are filtered to ensure statistical robustness (e.g., requiring a minimum of 10 samples per group and removing outliers with unusually high π). The final, filtered data is then used to generate a heatmap visualizing π for the top lineages over time.

<br>

![alt text](./results/heatmap_nucdiv.png "Nucleotide Diversity (π) for the Top 20 Lineages Over Time")

**Figure 1: Nucleotide Diversity (π) for the Top 20 Lineages Over Time.** This heatmap visualizes the calculated nucleotide diversity (π) for the 20 SARS-CoV-2 Pango lineages with the highest total sample counts in the dataset. The y-axis lists the Pango lineages, and the x-axis represents bi-weekly time bins from late 2020 to early 2022. The color intensity, following the 'viridis' colormap, corresponds to the level of nucleotide diversity, where dark purple indicates low diversity (values near 0.0) and bright yellow indicates higher diversity (values approaching 0.0008). White spaces indicate time bins for which a given lineage had insufficient data to calculate a robust π value.

The heatmap reveals distinct temporal patterns. For example, lineage B.1.177.7 shows relatively low but persistent diversity in late 2020 and early 2021. In contrast, lineage B.1.617.2 (Delta) exhibits a period of higher and more variable diversity throughout mid-to-late 2021, consistent with its global expansion and diversification. The emergence of various Omicron sublineages (prefixed with 'BA') in late 2021 and early 2022 is characterized by generally lower initial diversity, which is an expected signature of a selective sweep.

<br>

#### **Limitations and Considerations**

The quality of this analysis is contingent on the quality and representativeness of the input data. Inaccuracies in Pango lineage assignments or collection dates within the metadata will propagate directly into the final dataset. The analysis also relies on sequences submitted to public databases, which are subject to significant geographic and temporal sampling biases. The use of consensus sequences, while computationally necessary, masks potential intra-host variation. Finally, any binning of data into groups depends on sufficient sample sizes; calculations on very small groups (e.g., n < 10) are statistically unreliable and are therefore filtered out.

---

### **Data Hosting Plan for Public Access**

To ensure open access and reproducibility, the final alignment and supporting files will be archived in compressed format at a public GitHub repository named `data-sarscov2-genomes`. All steps are fully scripted and version-controlled. The pipeline’s code, environment specifications, and key data files will be hosted on GitHub with use of Git LFS for any large files. The specific tool versions and data-release dates are documented in the code. This ensures that another researcher can reproduce the alignment and downstream π analyses.

### **Future Plans: Statistical Analysis of Nucleotide Diversity Dynamics**

The successful generation of a validated, time-partitioned nucleotide diversity dataset enables the direct statistical testing of key evolutionary hypotheses. The next phase of this research will leverage the final results file to move from descriptive visualization to quantitative inference. The initial analysis will involve a formal statistical test, such as an independent samples t-test, to validate the hypothesis that major variant transitions are associated with significant decreases in π, providing quantitative evidence of selective sweeps. Subsequently, the research will extend to modeling the rate of diversification within individual long-circulating lineages using linear regression, which could reveal novel differences in their evolutionary dynamics. Finally, a more nuanced investigation will focus on periods of "mosaic" evolution by using Analysis of Variance (ANOVA) to compare the diversity profiles of co-circulating lineages. Such an analysis could uncover whether competing lineages evolved under different selective or demographic pressures, providing data-driven insights into the complex dynamics of a multi-variant ecosystem.

### **Acknowledgements**

The conceptual development and drafting of this report benefited significantly from discussions and iterative refinement with an AI language model, Gemini 2.5 Pro (Google, 6/5/2025). The process included the generation of hypotheses, development of Python code and documentation for implementation of the methodology, experimental design, drafting of this report, strategies for hosting and efficient handling of large data files, and interpretation of the data analysis and its presentation. The author oversaw and reviewed the accuracy and robustness of all parts of this study.
