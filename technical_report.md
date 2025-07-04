## **Technical Report: A Pipeline for Scalable Analysis of SARS-CoV-2 Nucleotide Diversity**

### **Abstract**

This report provides a comprehensive guide to the `nucdiv_v4.py` workflow, designed to calculate genome-wide nucleotide diversity (π) for SARS-CoV-2 using the UShER mutation-annotated phylogeny. Drawing upon the theoretical framework outlined in the accompanying research proposal, this document details each stage of the pipeline: data acquisition, dependency management, clade-based subsampling, and parallelizable consensus sequence generation using `matUtils` and `bcftools`. The result is a single, concatenated FASTA alignment, structured for direct use in downstream population genetics software. We present detailed code annotations, expected data outputs, and a robust data hosting strategy using a GitHub repository with Git LFS for any large files to ensure full reproducibility and community access.

---

### **Introduction**

#### **Background on this Report**

Understanding SARS-CoV-2 evolution requires linking genetic diversity to population processes. SARS-CoV-2 variants emerge through mutation, selection, recombination, and chance (genetic drift) under immune pressure. For example, while most viral mutations are neutral or deleterious, a few confer enhanced transmissibility or immune escape. These advantageous mutations (e.g., in Spike) can fix rapidly, driving **selective sweeps** that purge diversity. Conversely, when many lineages co-circulate (a “mosaic” regime), diversity remains high and recombination may shuffle genomes without an immediate sweep. Spatial spread also matters: when a variant spreads to new regions, founder bottlenecks tend to reduce diversity away from the origin.

To quantify these dynamics, we compute **nucleotide diversity (π)**—the average pairwise genetic differences per site among sampled genomes. Under neutral evolution, π equilibrates at approximately 2Nₑμ (where Nₑ is the effective population size and μ the per-site mutation rate). Thus, π serves as an empirical gauge of the viral population’s “effective size” and past bottlenecks.

In this report, we describe a scalable pipeline (`nucdiv_v4.py`) to build multi-FASTA alignments for π analysis, and we integrate population-genetic theory to interpret the outputs. The accompanying proposal establishes the rationale for a clade-based subsampling strategy, which balances representative sampling of major lineages against computational feasibility. By selecting up to 1,000 genomes from each of the 250 most prevalent Pangolin lineages, this approach mitigates biases from uneven sequencing efforts while retaining the statistical power needed to detect significant patterns in nucleotide diversity (π).

The core of this method is the UShER mutation-annotated tree, which provides both the global phylogenetic context and per-sample variant calls in a compressed protobuf format. The `nucdiv_v4.py` script leverages `matUtils` for rapid variant extraction and `bcftools` for consensus sequence generation. The final output is a single, concatenated FASTA file where each sequence header is annotated with its sample identifier and Pango lineage (e.g., `>sample_id|clade_name`), making it immediately suitable for downstream analysis with tools like `scikit-allel`, R, or for phylodynamic inference with BEAST.

In this report, Section 2 details each component of the Python script, explaining its function and rationale. Section 3 discusses expected runtimes and output formats. Finally, Section 4 outlines a version-controlled data hosting plan on GitHub to promote open science.

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

#### **Methods for Informatics Analysis and Interpretation**

The workflow proceeds from environment setup to final file cleanup. Each step is designed to be modular and reproducible.

**Data Flow Overview**

The pipeline follows a linear data flow, transforming raw public data into a final, analysis-ready alignment.

```
[Input Data] -> [Cell 1-3: Setup & Acquisition] -> [Cell 4: Subsampling] -> [Cell 5: Execution Engine] -> [Final Alignment]
  - .pb file       - Install tools               - Generate control file     - Extract VCFs           - all_clades_aligned.fas
  - .tsv metadata  - Download data               - Generate sample lists     - Build consensus
  - .fa reference                                                          - Append to FASTA
```

The workflow uses the public UShER mutation-annotated phylogeny (protobuf format) as input. Key steps (and code cells) include:

*   **Environment setup:** Install bioinformatics tools (UShER `matUtils`, `bcftools`) via Conda, and define file paths (UShER `.pb` file, metadata TSV, reference FASTA, etc.).
*   **Master sample extraction (Cell 1):** Extract all sample IDs from the UShER tree using `matUtils summary`, ensuring our analysis considers only genomes present in both the tree and metadata.
*   **Data acquisition (Cells 2–3):** Download a date-stamped UShER phylogeny, metadata, and reference genome for reproducibility. Using fixed release dates (not “latest”) ensures results can be exactly reproduced.
*   **Clade-based subsampling (Cell 4):** Group all samples by Pango lineage, select the top 250 lineages by prevalence, and randomly subsample up to 1,000 genomes per clade. This stratified design balances lineage representation while keeping data sizes manageable. Sample IDs are written to `samples_<clade>.txt` files, and a control file lists the clade names. Any clade names that do not represent valid clades, such as `Unassigned`, are removed.
*   **Consensus sequence generation (Cell 5):** For each subsampled clade, use `matUtils extract` to pull its variant calls into a VCF. We bgzip/index the VCF, then loop over sample IDs to run `bcftools consensus`, applying each genome’s variants to the reference to reconstruct its full sequence. The resulting FASTA records are appended to `all_clades_aligned.fas`, with headers like `>sample_id|CLADENAME` for easy downstream parsing.
*   **Verification and cleanup (Cell 6):** After all clades, the script checks that `all_clades_aligned.fas` exists and is non-empty, counts the total sequences, and removes intermediate files.

These steps produce a single concatenated FASTA alignment containing one consensus sequence per sample, annotated by lineage. The pipeline code is fully parameterized (lineage count, samples per clade, seed) and version-controlled, with data stored via Git LFS to ensure reproducibility of inputs, code, and outputs.

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

The commands shown here ensure a consistent analytical environment. `conda install` fetches the specified versions of UShER and bcftools. The `wget` commands download specific, dated releases of the UShER phylogeny and metadata. Using an archived version instead of "latest" is critical for reproducibility, as the live data files can be updated daily.

#### **Cell 4: Clade-Based Subsampling**

This cell implements the core sampling strategy defined in the proposal.

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
            # Determine sample size, ensuring not to exceed available samples
            num_to_sample = min(SAMPLES_PER_CLADE, len(group_df))
            subsampled_group = group_df.sample(n=num_to_sample, random_state=42)

            # Sanitize clade name for use as a valid filename
            safe_clade_name = "".join(c if c.isalnum() else "_" for c in clade_name)

            # Write sample list and update the control file
            subsampled_group['strain'].to_csv(f'samples_{safe_clade_name}.txt', index=False, header=False)
            f_control.write(f"{safe_clade_name}\n")
```

The script first defines the sampling parameters: the top 250 clades and up to 1,000 samples per clade. It loads the master sample list and filters the metadata to include only common entries. The `groupby()` and `nlargest()` methods efficiently identify the most prevalent clades. The script then iterates through these major clades. For each, it performs a reproducible random sample (`random_state=42`) of the specified size. The clade name is sanitized to remove special characters (e.g., `BA.2.75` becomes `BA_2_75`), which ensures it can be used safely in filenames. Finally, it writes the list of subsampled strain IDs to a dedicated `samples_*.txt` file and adds the sanitized clade name to the `clades_to_process.txt` control file.

#### **Cell 5: Consensus Sequence Generation Loop**

This shell script block is the execution engine of the pipeline, iterating through the work plan generated in Cell 4.

```bash
#!/bin/bash
set -e # Exit immediately if a command fails

# Define file paths
PB_FILE="public-latest.all.masked.pb"
REF_FASTA="NC_045512v2.fa"
CLADES_CONTROL_FILE="clades_to_process.txt"
FINAL_ALIGNED_FASTA="all_clades_aligned.fas"

# Initialize the final output file, but do not delete if it exists (for restarts)
touch "$FINAL_ALIGNED_FASTA"

# Read each sanitized clade name from the control file
while read -r safe_clade_name; do
    CLADE_SAMPLE_FILE="samples_${safe_clade_name}.txt"
    TEMP_VCF="temp_${safe_clade_name}.vcf"

    # 1. Extract variants for the current clade's samples into a VCF file
    matUtils extract -i "$PB_FILE" --samples "$CLADE_SAMPLE_FILE" -v "$TEMP_VCF"

    # 2. Compress and index the VCF for efficient access by bcftools
    bgzip --force "$TEMP_VCF"
    bcftools index --force --tbi "${TEMP_VCF}.gz"

    # 3. Loop through each sample ID to generate its consensus sequence
    while read -r sample_id; do
        bcftools consensus -f "$REF_FASTA" --sample "$sample_id" "${TEMP_VCF}.gz" \
          | sed "1s#.*#>${sample_id}|${safe_clade_name}#" >> "$FINAL_ALIGNED_FASTA"
    done < "$CLADE_SAMPLE_FILE"

    # 4. Clean up intermediate files for this clade to save space
    rm "${TEMP_VCF}.gz" "${TEMP_VCF}.gz.tbi"
done < "$CLADES_CONTROL_FILE"
```

The loop reads one clade name at a time from the control file. For each clade, it executes a four-step process:
1.  **`matUtils extract`**: Isolates the genetic variants specific to the samples in the current clade and writes them to a temporary VCF file.
2.  **`bgzip` & `bcftools index`**: Compresses and indexes the VCF. This is a mandatory prerequisite for `bcftools`, as it allows for rapid lookups of specific sample data.
3.  **`bcftools consensus`**: A nested loop reads each sample ID. For each one, `bcftools` reconstructs the full genome sequence by applying that sample's unique variants to the reference FASTA. The `sed` command modifies the FASTA header to include both the sample ID and its clade, then appends the result to the final output file.
4.  **`rm`**: Removes the temporary VCF files to conserve disk space.

#### **Cell 6: Verification and Cleanup**

The final cell verifies the output and cleans the working directory.

```bash
# Verify that the final alignment file was created and is not empty
if [ ! -s "$FINAL_ALIGNED_FASTA" ]; then
    echo "Error: Final alignment file not found or is empty." >&2
    exit 1
fi

# Count sequences to provide a summary of the run
SEQUENCE_COUNT=$(grep -c ">" "$FINAL_ALIGNED_FASTA")
echo "Verification complete. Final file contains $SEQUENCE_COUNT sequences."

# Remove all intermediate sample lists and the control file
rm -f clades_to_process.txt samples_*.txt
```

This script performs two checks: first, that the final FASTA file exists and has a size greater than zero (`-s` flag); second, it counts the number of headers (`>`) to report the total number of sequences generated. Finally, it removes all intermediate `samples_*.txt` files and the control file, leaving a clean directory with only the primary data and final result.

---

### **Part I: Consensus Sequence Generation from the UShER Tree**

This part details the `nucdiv_v4.py` workflow, designed for the initial, large-scale generation of a multiple sequence alignment from the raw UShER protobuf file.

#### **Methods for Interpretation of Nucleotide Diversity Patterns**

Once π is computed from the final alignment (e.g., for global and regional time-bins), we interpret the patterns against our theoretical expectations. Key anticipated signals and the methods to detect them include:

*   **Co-circulation/mosaic diversity:** Periods like late 2022, when multiple Omicron sublineages co-existed, should show **elevated and relatively stable π**. In these regimes, recombinants like XBB may appear, but their initial spread does not immediately suppress diversity. We expect the addition of XBB to **maintain** (or even slightly increase) π, rather than drop it. In practice, this means the pipeline alignment from, say, October–December 2022 should produce a high π estimate compared to preceding or succeeding intervals.
*   **Antigenic leaps:** The sudden emergence of a deep-branch variant (like Omicron BA.1 or BA.2.86) effectively adds a highly divergent genome to the sample. This can cause a *transient increase* in π as the gene pool’s average divergence jumps. We will check for such events by comparing π before and after known leap introductions; an unexplained diversity “spike” may flag a cryptic lineage’s appearance.
*   **Spatial founder effects:** When a sweep lineage spreads between regions, we expect a diversity gradient. For a recent variant, the origin region (with many introductions) should show higher π, whereas remote regions (seeded by few genomes) should show lower π early on. Over time, as the variant accumulates local mutations, π in sinks will rise. By computing π by region (e.g., origin vs. Europe vs. Asia) at each time-point, we can detect these gradients. Any systematic regional differences (tested via ANOVA or mixed models) would support founder-effect dynamics.
*   **Selective sweeps:** We expect sharp falls in π coinciding with the rise of a new dominant variant. For example, historical data show global diversity plunged when Alpha and Delta spread; similarly, the JN.1 sweep should cause a steep diversity collapse. In our pipeline output, a plot of π over time should reveal a pronounced minimum during JN.1’s emergence, reflecting the loss of polymorphism as one lineage took over.
*   **Partitioning total vs. lineage-specific π:** In mosaic phases, total diversity arises from deep splits between co-circulating lineages. We will explicitly compare the global π (`π_total`) in a window to the π computed *within* each major lineage (e.g., π_XBB, π_BQ.1). If multiple divergent lineages are present, we expect π_total ≫ π_lineage, because most diversity is between lineages. A much larger π_total indicates that the high diversity is driven by differences *across* lineages, a quantitative signature of a mosaic regime.
*   **BEAST phylodynamics (Optional):** For targeted time-frames or clades, BEAST is applicable for testing subsampled sequences. For instance, a BEAST skyline on JN.1 genomes could reveal an abrupt contraction in inferred Nₑ at sweep onset. Matching such Nₑ(t) drops to π crashes strengthens causal inference. Similarly, comparing Nₑ(t) for the ancestral XBB clade versus its major descendant (XBB.1.5) could test if later adaptations (rather than recombination alone) drove a sweep.

Taken together, these analyses allow us to translate pipeline outputs (the alignment and π values) into biological conclusions. For example, seeing a π collapse in the alignment-derived data for September 2023 would confirm the selective sweep of JN.1, while sustained high π through 2022 would illustrate the Omicron co-circulation era. By moving beyond qualitative variant narratives to quantitative π measurements, we validate the sweep-versus-mosaic model of SARS-CoV-2 evolution.

---

#### **Results of the Pipeline and its Methods**

##### **Expected Output Format**

The primary output is `all_clades_aligned.fas`, a multi-FASTA file. Each record follows a standardized format, with the header containing the original strain identifier and the sanitized Pango lineage, separated by a pipe (`|`).

```fasta
>USA/AZ-544/2020|BA_2_75
NNNNACGTACGT...
>another/sample/id|BQ_1_1
NNNNACGTACGT...
```

The aligned sequence data is further processed downstream for validity. The sequence names, sequence lengths, and their terminal regions are verified against their expected values and any data samples that are not verifiable by these criteria are omitted from the data file and  saved to a file for inspection.

##### **The Final Aligned Sequence Dataset and Its Analytical Utility**

The `all_clades_aligned.fas` file is the central product of this workflow. It is a concatenated multiple sequence alignment where each sequence is a full-length consensus genome, reconstructed by applying sample-specific variants from the UShER tree to a common reference. This whole-genome alignment is a powerful resource for several downstream applications:

1.  **Partitioned Evolutionary Analysis:** The SARS-CoV-2 genome contains different functional elements (e.g., Spike gene, ORF1ab, non-coding regions) that are subject to different evolutionary pressures. The whole-genome alignment allows for partitioning the data by gene or codon position. This enables the application of different evolutionary models to different genomic regions—for instance, using a codon-based model to detect positive selection in the Spike protein while applying a simpler nucleotide model to non-coding regions.

2.  **Phylogenetic Reconstruction:** This alignment serves as direct input for phylogenetic inference software such as IQ-TREE, RAxML, and BEAST. Reconstructing a high-resolution phylogeny from this dataset allows for detailed visualization of the evolutionary relationships among all sampled lineages. The resulting tree is fundamental for subsequent analyses.

3.  **Advanced Phylodynamic and Phylogeographic Inference:** A time-calibrated phylogeny derived from this alignment can be used to infer complex population dynamics. This includes estimating changes in effective population size over time (phylodynamics) and reconstructing the geographic spread of viral lineages (phylogeography). These methods provide quantitative insights into the processes driving the patterns of nucleotide diversity observed.

4.  **Ancestral State Reconstruction (ASR):** With a robust phylogeny, it is possible to perform ASR to infer the genetic sequences of ancestral viruses at key nodes in the tree. For example, one can reconstruct the genome of the most recent common ancestor of the JN.1 lineage. This allows for the precise identification of mutations that occurred on the evolutionary path leading to a new variant of concern, thereby providing direct evidence of the specific genetic changes that may have driven its emergence and success.

The pipeline was also designed for the limits of computation that is available to a typical researcher. Therefore, the computational power and memory requirements are for a workstation level system. In particular, the above aligned genomic sequence data requires approximately 7 gigabytes of disk space. Therefore, the methods and code of this study are adapted for higher efficiency, including the methods of file transfer, its storage on disk, and use of a memory-mapped binary file format for the input and output operations of the sequence data.

This entailed a modular-like approach in the steps for data collection, processing, and analysis. This was a practical concern in development of the informatics pipeline to perform these tasks since additions to it can lead to an an ever-increasing level of code complexity both within and between tasks.

##### **Runtime and Scalability**

Based on the proposal's benchmarking, processing one clade of 1,000 samples takes approximately three minutes on a standard compute node. For 250 clades, the total runtime is estimated to be over 12 hours. The workflow is parallelizable; the `clades_to_process.txt` file can be split into multiple chunks and run concurrently on different machines to reduce the wall-clock time.

##### **Limitations and Considerations**

The quality of this analysis is contingent on the quality and representativeness of the input data. Inaccuracies in Pango lineage assignments or collection dates within the metadata will propagate directly into the final dataset. The analysis also relies on sequences submitted to public databases, which are subject to significant geographic and temporal sampling biases. As sequencing efforts are not uniform globally, the resulting dataset is not a truly random sample of the viral population and may underrepresent diversity in certain regions or time periods.

Furthermore, the use of consensus sequences, while computationally necessary for an analysis of this scale, masks potential intra-host variation. The diversity metrics calculated therefore reflect variation *between* circulating viral populations, not the diversity present within a single infected individual. The consensus sequences are also generated directly from the UShER-provided VCFs and reference genome; they are not re-aligned or masked for problematic regions (e.g., hypervariable sites or sequencing artifacts), which may be a necessary preprocessing step for certain downstream analyses.

Finally, any binning of data into groups depends on sufficient sample sizes. The calculation of diversity metrics on very small sample sizes (e.g., n < 10) can produce statistically unreliable results. Standard practice is therefore applied to filter out these groups before performing the main calculations and visualizations.

---

### **Part II: Analysis of Nucleotide Diversity from Aligned Sequences**

The `nucdiv_stats.py` script, detailed in this section, processes the `all_clades_aligned.fas` file generated in Part I. Its primary goal is to calculate nucleotide diversity (π) across different time periods and viral lineages. This script is designed to handle large (multi-gigabyte) datasets efficiently on standard computer hardware by avoiding loading entire files into memory simultaneously. The workflow within `nucdiv_stats.py` (represented by the "Part 1" and "Part 2" Python script segments provided by the user) can be broken down into several key stages:

1.  **Initial Data Acquisition and Verification (Script Part 1)**:
    *   Combines FASTA file parts (if the alignment was split) from a persistent storage like Google Drive.
    *   Verifies each sequence in the combined FASTA file for correct header format (expecting `sample_id|clade_name|...|...` format, specifically checking for 3 pipe delimiters) and consistent sequence length.
    *   Outputs three critical files:
        *   `verified_sequences.fas`: A clean FASTA file containing only valid sequences.
        *   `verification_log.csv`: A log of discarded sequences and reasons for discarding.
        *   `verified_sequence_order.txt`: A text file listing the sequence IDs (just the `sample_id` part of the header) in the exact order they appear in `verified_sequences.fas`. This index is crucial for later metadata alignment.
2.  **Pre-processing for Efficient Analysis (Script Part 2, Cell 2)**:
    *   Converts `verified_sequences.fas` into a NumPy memory-mapped binary file (`alignment.mmap`). This allows for efficient, random access to sequence data without loading the entire alignment into RAM.
3.  **Memory-Efficient Metadata Alignment (Script Part 2, Cell 3)**:
    *   Aligns a comprehensive metadata file (e.g., `public-latest.metadata.tsv`) with the verified sequence data.
    *   Reads the large metadata file in chunks.
    *   Filters each chunk to retain only entries corresponding to sequence IDs found in `verified_sequence_order.txt`.
    *   Concatenates filtered chunks, de-duplicates entries by 'strain' ID (keeping the first occurrence).
    *   Re-indexes the resulting metadata DataFrame to precisely match the order in `verified_sequence_order.txt`.
    *   Creates temporal partitions by converting dates to bi-weekly time bins (`time_bin`).
4.  **Partitioned Nucleotide Diversity (π) Calculation (Script Part 2, Cell 4)**:
    *   Groups the aligned metadata by `pangolin_lineage` and `time_bin`.
    *   For each group (if it contains at least 2 samples):
        *   Retrieves the corresponding sequences from the `alignment.mmap` file using the indices from the grouped metadata.
        *   Calculates π for this sequence slice. The π calculation involves summing pairwise differences at each site and normalizing by the number of pairs and sequence length.
    *   Stores results (lineage, time_bin, π, number of samples) in a list.
5.  **Results Output, Diagnostic, Filtering, and Visualization (Script Part 2, subsequent cells)**:
    *   Saves the raw π calculation results to `nucleotide_diversity_results.csv`.
    *   Includes an optional diagnostic step to investigate groups with unusually high π values (e.g., > 0.001), identifying the most divergent sequence pair within such groups.
    *   Filters the results to ensure statistical robustness (e.g., requiring a minimum number of samples like 10, and removing groups with π above the defined threshold). Saves this as `nucleotide_diversity_results_filtered.csv`.
    *   Generates a heatmap using `seaborn` and `matplotlib` to visualize π for the top N lineages over time, providing a comprehensive overview of diversity dynamics.

The script segments provided illustrate a robust pipeline for moving from a large, raw sequence alignment to quantitative evolutionary insights.

#### **Detailed Script Breakdown and Methods**

Below, sections from the `nucdiv_stats.py` script (as provided by the user in Jupyter-like cells) are integrated and explained.

**Part 1: Standalone Verification and Filtering Script**

This initial part of the workflow focuses on preparing and validating the input FASTA file.

*   **File Handling and Setup**: The script begins by setting up access to Google Drive (if used for storage) and defining paths for input and output files. It can assemble a complete FASTA file from multiple parts if the alignment was previously split. The following Python snippet illustrates local file configuration. Shell commands would handle file assembly from parts (e.g., `gunzip -c *.fas.gz > all_clades_aligned_combined.fas` after copying parts).

    ```python
    # Part 1: Configuration (Illustrative Python part)
    import os
    # import csv # csv is used in functions below
    # from google.colab import drive # If using Colab

    # # Mount Google Drive (if applicable)
    # # drive.mount('/content/drive')
    # # DRIVE_PATH="/content/drive/MyDrive/"

    # --- Configuration for local execution ---
    # Assuming combined FASTA is already created by shell commands if it was in parts
    INPUT_FASTA = "all_clades_aligned_combined.fas" # Combined raw FASTA
    CLEAN_FASTA_OUTPUT = "verified_sequences.fas"
    LOG_FILE_OUTPUT = "verification_log.csv"
    ORDER_INDEX_OUTPUT = "verified_sequence_order.txt" # Critical for metadata alignment
    ```

*   **FASTA Iterator and Verification Logic**: The `fasta_iterator` function efficiently reads large FASTA files. The `verify_and_filter_fasta` function applies validation rules: it checks for correctly formatted headers (expecting 3 pipe characters, implying 4 fields like `>strain|lineage|date|country`) and consistent sequence lengths. Valid sequences are written to `verified_sequences.fas`, discards are logged in `verification_log.csv`, and crucially, the primary sequence IDs (the part of the header before the first pipe) of valid sequences are written in order to `verified_sequence_order.txt`.

    ```python
    import csv # Make sure csv is imported

    def fasta_iterator(file_path):
        """A generator to read a FASTA file one sequence at a time."""
        with open(file_path, 'r') as f:
            header, sequence_lines = None, []
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith('>'):
                    if header: yield header, "".join(sequence_lines)
                    header = line[1:] # Store header without '>'
                    sequence_lines = []
                else:
                    sequence_lines.append(line)
            if header: yield header, "".join(sequence_lines)

    def verify_and_filter_fasta(input_path, clean_path, log_path, order_path):
        """
        Streams a FASTA file, verifies each sequence, writes valid ones to a new
        file, logs discards, and creates an ordered index of valid IDs.
        """
        print(f"Starting verification of '{input_path}'...")
        if not os.path.isfile(input_path):
            print(f"FATAL: Input file '{input_path}' not found.")
            return 0, 0 # Return counts of processed and written

        processed, written, errors = 0, 0, 0
        ref_len = None # Reference length is determined by the first valid sequence

        with open(clean_path, 'w') as f_clean, \
             open(log_path, 'w', newline='') as f_log, \
             open(order_path, 'w') as f_order:

            log_writer = csv.writer(f_log)
            log_writer.writerow(['full_header', 'error_reason'])

            for header, sequence in fasta_iterator(input_path):
                processed += 1
                error_reason = None
                # Extract the primary sequence ID (part before the first pipe if pipes exist)
                seq_id_for_order = header.split('|')[0] if '|' in header else header

                if header.count('|') != 3: # Expecting format like >id|clade|date|country
                    error_reason = "Malformed header (expected 3 pipe characters)"
                else:
                    if ref_len is None:
                        ref_len = len(sequence)
                        if ref_len == 0:
                            error_reason = "Zero length sequence"
                            ref_len = None # Reset so next sequence can try to set it
                    elif len(sequence) != ref_len:
                        error_reason = f"Incorrect length (is {len(sequence)}, expected {ref_len})"

                if error_reason:
                    errors += 1
                    log_writer.writerow([header, error_reason])
                else:
                    written += 1
                    f_clean.write(f">{header}\n{sequence}\n")
                    f_order.write(f"{seq_id_for_order}\n")

        print("\n--- Verification and Filtering Complete ---")
        print(f"Total sequences processed: {processed}")
        print(f"Sequences written to new file: {written}")
        print(f"Sequences discarded (errors): {errors}")
        # (Other print statements for file paths)
        return processed, written

    # --- Main Execution (Illustrative Call) ---
    # verify_and_filter_fasta(INPUT_FASTA, CLEAN_FASTA_OUTPUT, LOG_FILE_OUTPUT, ORDER_INDEX_OUTPUT)
    ```
    This verification step is critical. The `verified_sequence_order.txt` file ensures correct alignment between sequence data and metadata in subsequent steps.

**Part 2: Partitioned Nucleotide Diversity Script**

This part takes the verified data from Part 1 and calculates π.

*   **Setup and Data Acquisition (Cell 1 of Script Part 2)**:
    Acquires necessary files: `verified_sequences.fas` (gzipped), `verified_sequence_order.txt`, and the main `public-latest.metadata.tsv`. Files are typically copied from persistent storage (like Google Drive) to the local Colab environment and decompressed if needed.

    ```python
    # Part 2, Cell 1: Setup and Data Acquisition (Illustrative)
    import numpy as np
    import pandas as pd
    # import os # Assuming os is imported

    # --- Configuration ---
    # DRIVE_PATH = "/content/drive/MyDrive/"
    LOCAL_DATA_PATH = "./"
    CLEAN_FASTA_GZ_FROM_DRIVE = "verified_sequences.fas.gz" # Output of Part 1
    ORDER_INDEX_FILE_FROM_DRIVE = "verified_sequence_order.txt" # Output of Part 1
    METADATA_TSV_FROM_DRIVE = "public-latest.metadata.tsv"

    CLEAN_FASTA_LOCAL = os.path.join(LOCAL_DATA_PATH, "verified_sequences.fas")
    ORDER_INDEX_LOCAL = os.path.join(LOCAL_DATA_PATH, "verified_sequence_order.txt")
    METADATA_TSV_LOCAL = os.path.join(LOCAL_DATA_PATH, "public-latest.metadata.tsv")
    NUMPY_BINARY_FILE = os.path.join(LOCAL_DATA_PATH, "alignment.mmap")

    # --- File Acquisition (Illustrative shell commands if using Colab) ---
    # print("Acquiring input files from Google Drive...")
    # !cp "{DRIVE_PATH}{CLEAN_FASTA_GZ_FROM_DRIVE}" "{LOCAL_DATA_PATH}"
    # !cp "{DRIVE_PATH}{ORDER_INDEX_FILE_FROM_DRIVE}" "{LOCAL_DATA_PATH}"
    # !cp "{DRIVE_PATH}{METADATA_TSV_FROM_DRIVE}" "{LOCAL_DATA_PATH}"
    # !gunzip -f "{LOCAL_DATA_PATH}{CLEAN_FASTA_GZ_FROM_DRIVE}" # Results in CLEAN_FASTA_LOCAL
    # print("All data ready for processing.")
    ```

*   **Pre-processing: FASTA to NumPy Binary (Cell 2 of Script Part 2)**:
    The `fasta_to_numpy_binary` function converts the text-based `verified_sequences.fas` into a numerical, memory-mapped NumPy array (`alignment.mmap`). Nucleotides (A,C,G,T,N) are mapped to integers (0,1,2,3,4). This binary format allows for much faster data access during π calculation.

    ```python
    # Part 2, Cell 2: Pre-processing (FASTA to Binary)
    # (Requires numpy as np, os)
    def fasta_to_numpy_binary(fasta_path, output_path):
        """Converts a FASTA alignment to a memory-mapped NumPy binary file."""
        print("Starting conversion from FASTA to memory-mapped binary format...")
        num_sequences, seq_length = 0, 0
        with open(fasta_path, 'r') as f: # First pass to get dimensions
            for line in f:
                if line.startswith('>'): num_sequences += 1
                elif not seq_length and num_sequences > 0: seq_length = len(line.strip())

        if num_sequences == 0 or seq_length == 0:
            print(f"Error: Found {num_sequences} sequences or seq_length {seq_length}. Cannot proceed.")
            return None, 0, 0
        print(f"Alignment dimensions: {num_sequences} sequences, {seq_length} sites.")

        mmap_array = np.memmap(output_path, dtype=np.uint8, mode='w+', shape=(num_sequences, seq_length))
        nt_map = str.maketrans("ACGTN", "01234") # N and any other char will map to 4 via logic below

        current_seq_index = 0
        with open(fasta_path, 'r') as f:
            for line in f:
                line_stripped = line.strip().upper()
                if not line_stripped or line_stripped.startswith('>'): continue

                translated_sequence = line_stripped.translate(nt_map)
                # Ensure all chars are 0-4; map others to 4 (N)
                numeric_sequence = np.array([int(c) if c in '01234' else 4 for c in translated_sequence], dtype=np.uint8)

                if len(numeric_sequence) == seq_length:
                    mmap_array[current_seq_index, :] = numeric_sequence
                    current_seq_index += 1
                # else: error, should have been caught by verify_and_filter_fasta

        mmap_array.flush()
        # It's good practice to delete the memmap object to ensure the file is closed,
        # especially if it will be reopened later (e.g., in 'r' mode).
        del mmap_array
        print(f"Successfully created binary file at: {output_path}")
        return output_path, num_sequences, seq_length

    # --- Main Execution (Illustrative Call) ---
    # binary_file_path, n_total_seq, n_total_len = fasta_to_numpy_binary(CLEAN_FASTA_LOCAL, NUMPY_BINARY_FILE)
    ```

*   **Memory-Efficient Metadata Alignment (Cell 3 of Script Part 2)**:
    This step aligns the potentially very large metadata TSV file with the sequence data. It reads the metadata in chunks (`CHUNK_SIZE`), filters for relevant strains using the IDs from `verified_sequence_order.txt`, concatenates, de-duplicates (keeping the first entry per strain), and then re-indexes the resulting DataFrame to match the exact order of sequences in the binary alignment. Finally, it creates a `time_bin` column (e.g., bi-weekly periods) from sample collection dates.

    ```python
    # Part 2, Cell 3: Metadata Alignment
    # (Requires pandas as pd, gc)
    # ORDER_INDEX_FILE = ORDER_INDEX_LOCAL # Path to verified_sequence_order.txt
    # METADATA_TSV = METADATA_TSV_LOCAL     # Path to public-latest.metadata.tsv
    # CHUNK_SIZE = 500000

    print("--- Aligning Metadata with Sequence Data (Memory-Efficient) ---")
    ordered_ids_set_for_filtering = set(pd.read_csv(ORDER_INDEX_LOCAL, header=None, names=['strain'])['strain'])

    cols_to_use = ['strain', 'pangolin_lineage', 'date', 'country'] # Customize as needed
    col_dtypes = {'strain': 'object', 'pangolin_lineage': 'category', 'date': 'object', 'country': 'category'}

    list_of_filtered_chunks = []
    with pd.read_csv(METADATA_TSV_LOCAL, sep='\t', usecols=cols_to_use, dtype=col_dtypes, chunksize=CHUNK_SIZE, on_bad_lines='warn') as reader:
        for chunk_df in reader:
            # Normalize strain ID in metadata (e.g. EPI_ISL_XXXXX|2023-01-01 -> EPI_ISL_XXXXX)
            # This must match how IDs were stored in verified_sequence_order.txt
            chunk_df['strain_normalized_for_join'] = chunk_df['strain'].str.split('|').str[0]
            filtered_chunk = chunk_df[chunk_df['strain_normalized_for_join'].isin(ordered_ids_set_for_filtering)]
            list_of_filtered_chunks.append(filtered_chunk)

    if not list_of_filtered_chunks:
        print("FATAL: No matching records found in metadata for sequence IDs.")
        df_meta_aligned = pd.DataFrame() # Assign empty df
    else:
        df_meta_concat_filtered = pd.concat(list_of_filtered_chunks, ignore_index=True)
        # del list_of_filtered_chunks; gc.collect()

        df_meta_deduplicated = df_meta_concat_filtered.drop_duplicates(subset=['strain_normalized_for_join'], keep='first')
        # del df_meta_concat_filtered; gc.collect()

        ordered_ids_list_for_reindex = pd.read_csv(ORDER_INDEX_LOCAL, header=None, names=['strain_normalized_for_join'])['strain_normalized_for_join'].tolist()
        df_meta_aligned = df_meta_deduplicated.set_index('strain_normalized_for_join').loc[ordered_ids_list_for_reindex].reset_index()
        # del df_meta_deduplicated, ordered_ids_list_for_reindex; gc.collect()

        df_meta_aligned['date'] = pd.to_datetime(df_meta_aligned['date'], errors='coerce')
        df_meta_aligned.dropna(subset=['date'], inplace=True)
        df_meta_aligned['time_bin'] = df_meta_aligned['date'].dt.to_period('2W') # Bi-weekly time bins

        # grouped_for_pi_calc = df_meta_aligned.groupby(['pangolin_lineage', 'time_bin'])
        # print(f"Metadata aligned and grouped. Shape: {df_meta_aligned.shape}, Groups: {len(grouped_for_pi_calc)}")
    ```

*   **Partitioned Nucleotide Diversity (π) Calculation (Cell 4 of Script Part 2)**:
    This cell iterates through each defined group (e.g., by `pangolin_lineage` and `time_bin`). For each group with at least two samples, it retrieves the corresponding sequence indices from `df_meta_aligned`, uses these to extract the sequence slice from `alignment.mmap`, and then calculates π using the `calculate_pi_for_slice` function. The function, as provided by the user, calculates π by summing, for each site, the proportion of differing pairs of nucleotides (considering only A,C,G,T) and averaging this over sites where comparison is possible.

    ```python
    # Part 2, Cell 4: Partitioned Pi Calculation
    # (Requires numpy as np, pandas as pd)
    # binary_file_path, n_total_seq, n_total_len should be available from fasta_to_numpy_binary
    # df_meta_aligned and grouped_for_pi_calc should be available from metadata alignment

    def calculate_pi_for_slice(alignment_slice_np):
        """Calculates pi for a given numpy array slice of an alignment."""
        num_seq_in_slice, len_seq = alignment_slice_np.shape
        if num_seq_in_slice < 2: return 0.0, 0

        total_pi_sum_across_sites = 0.0
        num_sites_compared = 0

        for k_site in range(len_seq): # Iterate over each site (column)
            column_at_site_k = alignment_slice_np[:, k_site]
            # Filter out 'N's or other non-ACGT characters (mapped to value 4)
            valid_bases_at_site_k = column_at_site_k[column_at_site_k < 4] # ACGT are 0,1,2,3

            n_k = len(valid_bases_at_site_k) # Number of valid sequences at this site
            if n_k < 2: continue # Need at least 2 sequences to compare at this site

            num_sites_compared += 1
            # Counts of A, C, G, T at this site
            counts_k = np.bincount(valid_bases_at_site_k, minlength=4) # minlength for A,C,G,T

            # Calculate pi for site k using Tajima's estimator: (n/(n-1)) * (1 - sum(p_i^2))
            sum_freq_sq_k = np.sum((counts_k / n_k)**2)
            pi_site_k = (n_k / (n_k - 1)) * (1 - sum_freq_sq_k)
            total_pi_sum_across_sites += pi_site_k

        # Average pi over sites where comparison was possible
        final_pi = total_pi_sum_across_sites / num_sites_compared if num_sites_compared > 0 else 0.0
        return final_pi, num_seq_in_slice # Return pi and number of samples used in this slice

    pi_results_data = []
    # if 'grouped_for_pi_calc' in locals() and binary_file_path and os.path.exists(binary_file_path):
    #     alignment_mmap_readonly = np.memmap(binary_file_path, dtype=np.uint8, mode='r', shape=(n_total_seq, n_total_len))
    #     for group_name_tuple, group_meta_df in grouped_for_pi_calc:
    #         if len(group_meta_df) < 2: continue # Skip if less than 2 samples

    #         indices_for_slice = group_meta_df.index.tolist() # Get original indices from df_meta_aligned
    #         current_seq_slice = alignment_mmap_readonly[indices_for_slice, :]

    #         pi_val, n_samps = calculate_pi_for_slice(current_seq_slice)
    #         # print(f"L: {group_name_tuple[0]}, T: {str(group_name_tuple[1])}, Pi: {pi_val:.6f}, N: {n_samps}")
    #         pi_results_data.append({
    #             'lineage': group_name_tuple[0],
    #             'time_bin': str(group_name_tuple[1]),
    #             'pi': pi_val,
    #             'n_samples': n_samps
    #         })
    #     del alignment_mmap_readonly # Close memmap
    #     results_df_raw = pd.DataFrame(pi_results_data)
    #     # results_df_raw.to_csv("nucleotide_diversity_results.csv", index=False)
    # else:
    #     print("Skipping Pi calculation: required data not found.")
    #     results_df_raw = pd.DataFrame()
    ```

*   **Diagnostic, Filtering, and Visualization (Subsequent Cells in Script Part 2)**:
    *   **Diagnostic Cell for High π Values**: Investigates groups with π > `PI_THRESHOLD` (e.g., 0.001). For each such group, it re-fetches the sequence slice and identifies the pair of sequences with the maximum nucleotide differences, printing their IDs. This helps pinpoint sources of unusually high diversity.
    *   **Filter Final Results**: The raw `nucleotide_diversity_results.csv` is filtered. Groups with `n_samples` < `MIN_SAMPLES_THRESHOLD` (e.g., 10) or `pi` > `PI_THRESHOLD` are removed. The cleaned data is saved as `nucleotide_diversity_results_filtered.csv`.
    *   **Generate Heatmap**: Uses `seaborn` and `matplotlib` to create a heatmap from `nucleotide_diversity_results_filtered.csv`. It typically shows the top `NUM_LINEAGES_TO_DISPLAY` (e.g., 20) by sample count on the y-axis, `time_bin` on the x-axis, and color intensity representing π. This provides a visual summary of diversity trends.

    ```python
    # --- (Illustrative) Diagnostic, Filtering, and Visualization ---
    # PI_THRESHOLD_DIAG = 0.001
    # MIN_SAMPLES_FILTER = 10
    # NUM_TOP_LINEAGES_HEATMAP = 20

    # # ... (Code for diagnostic cell using results_df_raw, df_meta_aligned, alignment_mmap_readonly)

    # if not results_df_raw.empty:
    #     results_df_filtered_final = results_df_raw[
    #         (results_df_raw['n_samples'] >= MIN_SAMPLES_FILTER) &
    #         (results_df_raw['pi'] <= PI_THRESHOLD_DIAG) # Using same threshold for consistency
    #     ].copy()
    #     # results_df_filtered_final.to_csv("nucleotide_diversity_results_filtered.csv", index=False)

    #     # ... (Code for heatmap generation using results_df_filtered_final)
    # else:
    #     print("Skipping final filtering and heatmap as no raw results were generated.")
    ```

This comprehensive pipeline allows for scalable and robust calculation of nucleotide diversity from large SARS-CoV-2 genomic datasets, enabling detailed investigation of viral evolutionary dynamics. The use of memory-mapping, chunked processing, and careful indexing are key to its efficiency.

<br>

![alt text](./results/heatmap_nucdiv.png "Nucleotide Diversity (π) for the Top 20 Lineages Over Time")

**Figure 1: Nucleotide Diversity (π) for the Top 20 Lineages Over Time.** This heatmap visualizes the calculated nucleotide diversity (π) for the 20 SARS-CoV-2 Pango lineages with the highest total sample counts in the dataset. The y-axis lists the Pango lineages, and the x-axis represents bi-weekly time bins from late 2020 to early 2022. The color intensity, following the 'viridis' colormap, corresponds to the level of nucleotide diversity, where dark purple indicates low diversity (values near 0.0) and bright yellow indicates higher diversity (values approaching 0.0008). White spaces indicate time bins for which a given lineage had insufficient data to calculate a robust π value.

The heatmap reveals distinct temporal patterns. For example, lineage B.1.177.7 shows relatively low but persistent diversity in late 2020 and early 2021. In contrast, lineage B.1.617.2 (Delta) exhibits a period of higher and more variable diversity throughout mid-to-late 2021, consistent with its global expansion and diversification. The emergence of various Omicron sublineages (prefixed with 'BA') in late 2021 and early 2022 is characterized by generally lower initial diversity, which is an expected signature of a selective sweep.
<br>

### **Data Hosting Plan for Public Access**

To ensure open access and reproducibility, the final alignment and supporting files will be archived in compressed format at a public GitHub repository named `data-sarscov2-genomes`.

All steps above are fully scripted and version-controlled. The pipeline’s code, environment specifications, and key data files (final FASTA) and URL links (UShER `.pb` tree, reference genome) will be hosted on GitHub with use of Git LFS for any large files. The specific tool versions (e.g., UShER 0.5.8, bcftools 1.18) and the exact data-release dates are documented in the code. This ensures that another researcher can reproduce the alignment and downstream π analyses. We also log the random seed for subsampling to allow identical sampling when needed.

By combining classical population-genetic reasoning with scalable genomic workflows, this pipeline provides both the **means** (alignment data) and the **interpretation** (theory of π under selection, drift, recombination) to rigorously analyze SARS-CoV-2 diversity. Expected patterns (diversity collapse at sweeps, elevated diversity during co-circulation, spatial gradients from founder effects) are grounded in theory, and any use of BEAST phylodynamic models may contextualize them by inferring the hidden Nₑ(t) trajectories (as suggested by Suchard et al., 2018). Ultimately, a sudden drop in π in surveillance data could itself serve as an early warning of an ongoing selective sweep, while persistently high π could signal an era of active recombination. This integrated approach—quantitative π calculation informed by evolutionary theory—yields a comprehensive picture of SARS-CoV-2’s population dynamics.

### **Future Plans: Statistical Analysis of Nucleotide Diversity Dynamics**

The successful generation of a validated, time-partitioned nucleotide diversity dataset enables the direct statistical testing of key evolutionary hypotheses. The next phase of this research will leverage the final results file to move from descriptive visualization to quantitative inference. The initial analysis will involve a formal statistical test, such as an independent samples t-test, to validate the hypothesis that major variant transitions are associated with significant decreases in π, providing quantitative evidence of selective sweeps. Subsequently, the research will extend to modeling the rate of diversification within individual long-circulating lineages using linear regression, which could reveal novel differences in their evolutionary dynamics. Finally, a more nuanced investigation will focus on periods of "mosaic" evolution by using Analysis of Variance (ANOVA) to compare the diversity profiles of co-circulating lineages. Such an analysis could uncover whether competing lineages evolved under different selective or demographic pressures, providing data-driven insights into the complex dynamics of a multi-variant ecosystem.

### **Acknowledgements**

The conceptual development and drafting of this report benefited significantly from discussions and iterative refinement with an AI language model, Gemini 2.5 Pro, (Google, 6/5/2025). The process included the generation of hypotheses, development of Python code and documentation for implementation of the methodology, experimental design, drafting of this report, strategies for hosting and efficient handling of large data files, and interpretation of the data analysis and its presentation. The author oversaw and and reviewed the accuracy and robustness of all parts of this study.
