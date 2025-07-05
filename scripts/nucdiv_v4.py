'''
**A Data-Flow Perspective**

*   **Cells 1-3: Raw Material Acquisition**
    *   **Input:** Public URLs.
    *   **Process:** Download and decompress data.
    *   **Output:** The three primary files (`public-latest.all.masked.pb`, `public-latest.metadata.tsv`,
          `NC_045512v2.fa`). These are the foundational assets for the entire analysis.

*   **Cell 4: The Planning and Organization Stage**
    *   **Input:** The raw `metadata.tsv` and a list of all samples from the `.pb` file.
    *   **Process:** This cell acts as a project manager. It reads the list of all available samples and the metadata.
          It then cross-references them to determine which samples belong to which clade.
    *   **Output:** It does not produce the final sequence data. Instead, it produces a *work plan* for the next cell.
          This plan consists of two parts:
        1.  `clades_to_process.txt`: A simple list of tasks (the clades to be processed).
        2.  `samples_*.txt` files: A series of small files, where each file contains the specific list of samples
            needed for one task.
    *   *The key dependency is that Cell 5 is entirely reliant on this work plan. Without it, Cell 5 would not know
         which clades to process or which samples belong to each clade.*

*   **Cell 5: The Execution Engine**
    *   **Input:** The work plan from Cell 4 (`clades_to_process.txt` and `samples_*.txt`) and the raw `.pb` file.
    *   **Process:** It reads the task list. For each task (each clade), it uses the corresponding sample list to
          extract the relevant data from the main `.pb` file and generate the consensus sequences.
    *   **Output:** It appends each newly generated sequence to a single, cumulative file: `all_clades_aligned.fas`.

*   **Cell 6: Final Verification and Cleanup**
    *   **Input:** The final `all_clades_aligned.fas` file and the intermediate work plan files.
    *   **Process:** It inspects the final output to ensure it was created and contains data. It then removes the
          now-unnecessary work plan files.
    *   **Output:** A clean directory containing only the final results.
'''

'''
### **1. Analysis of a Single Task (Clade "A")**

From the log, analysis of the time taken for the first clade:
*   **Number of Samples:** 1,257
*   **Step 1 (`matUtils extract`):** This took approximately 63 seconds (~1 minute). This step's duration depends
      on the size of the clade.
*   **Step 3 (`bcftools consensus` loop):** This is the most time-consuming part. The script executes this command
    **1,257 times** for this clade alone. While each individual execution is fast, the cumulative time is significant.

### **2. Estimation for the Entire Workflow**

The key insight for a large-scale estimation is that the total runtime will be majorly dependent on the `bcftools consensus`
loop, as it runs once for *every single sample*.

*   **Total Samples:** Your log from Cell 4 indicates you are processing **~8.4 million** common samples
      (`Found 8376055 common samples...`).
*   **Total `bcftools` Executions:** The inner loop will therefore run approximately 8.4 million times over the
      entire course of the job.

Let us make a conservative assumption that each `bcftools consensus` call takes, on average, 0.2 seconds.

*   **Estimated Time:** 8,400,000 samples × 0.2 seconds/sample = 1,680,000 seconds
*   **Conversion to Hours:** 1,680,000 seconds / 3,600 seconds/hour ≈ 467 hours
*   **Conversion to Days:** 467 hours / 24 hours/day ≈ **19.5 days**

This estimation does not include the cumulative time for the `matUtils extract` step, which would add many
more hours to the total.

### **Critical Consideration: Colab Session Limits**

Google Colab sessions have a maximum lifetime, which is typically around 12 hours for the free tier.
They are not designed for multi-day, uninterrupted computations. The connection will be terminated long
before this script can complete.
'''

'''
CLADE SAMPLING STRATEGIES

**Partitioning:** The script uses the `pandas.groupby('pangolin_lineage')` function.
This function partitions the entire list of samples into mutually exclusive groups.
The group for `BA.2` will contain only those samples explicitly labeled as `BA.2` in the metadata,
and the group for `BA.2.75` will contain only those samples labeled `BA.2.75`.

**Output:** This script then writes these distinct, non-overlapping lists of sample IDs into
separate files (e.g., `samples_BA_2.txt`, `samples_BA_2_75.txt`).

The `exclusive_count` in the `clades.txt` file is a pre-calculated summary of the size of these unique groups.
The algorithm uses this count to decide which clades to process, and the script's `groupby` logic ensures that
the resulting sample lists for those clades are indeed disjoint.

*   **For Nucleotide Diversity (π):** The primary goal is to obtain a stable and precise estimate of π.
In population genetics, sample sizes in the range of 50–100 are often considered sufficient to characterize
the diversity of a single population at a single point in time, as the precision of the estimate improves with
diminishing returns for larger sample sizes. The script's parameter of `SAMPLES_PER_CLADE = 100` aligns with
this convention.
'''

### **Cell 1: Environment Setup and Configuration (Python)**

# This cell prepares the environment, mounts Google Drive, and defines all file paths.

# ==============================================================================
# Cell 1: Environment Setup and Configuration
# ==============================================================================
# Installation of Conda for Google Colab
!pip install -q condacolab
import condacolab
condacolab.install()

# --- Core Imports and Drive Mounting ---
import os
import pandas as pd
import sys

print(f"Current working directory set to: {os.getcwd()}")

# --- Path Definitions ---
# Input Data Paths
PB_FILE = "public-latest.all.masked.pb"
METADATA_FILE = "public-latest.metadata.tsv"
REF_FASTA = "NC_045512v2.fa"
SAMPLES_TSV = "samples.txt" # A complete list of samples in the PB file

# Intermediate and Final Output Paths
CLADES_CONTROL_FILE = "clades_to_process.txt"
FINAL_ALIGNED_FASTA = "all_clades_aligned.fas"

# This command generates the master sample list from the PB file.
# It is a prerequisite for the Python logic in Cell 4.
print("Generating master sample list from the protobuf file...")
!matUtils summary -i {PB_FILE} -s {SAMPLES_TSV}
print("Master sample list created.")

# The following would be uncommented in Colab, but in its own Cell
#  of shell code only (not mixed with Python code).
# Commented out IPython magic to ensure Python compatibility.
# ### **Cell 2: Install Bioinformatics Dependencies (Shell)**
# 
# # This cell uses Conda to install the required bioinformatics tools.
# 
# %%shell
# # ==============================================================================
# # Cell 2: Install Bioinformatics Dependencies
# # ==============================================================================
# echo "Installing bioinformatics tools via Conda..."
# conda install -y -c bioconda -c conda-forge usher bcftools
# echo "Installation complete."

# Commented out IPython magic to ensure Python compatibility.
# %%shell
# # ==============================================================================
# # Cell 3: Acquire a Consistent, Archived Data Set
# # ==============================================================================
# # Using a specific dated release (e.g., from June 29, 2025) to ensure data consistency.
# # The "latest" links can be out of sync.
# RELEASE_DATE="2025/06/29"
# BASE_URL="http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/${RELEASE_DATE}/"
# PB_FILE_NAME="public-${RELEASE_DATE//\//-}.all.masked.pb"
# METADATA_FILE_NAME="public-${RELEASE_DATE//\//-}.metadata.tsv"
# # CLADE_FILE_NAME="cladeToPublicName.tsv"
# # LINEAGE_FILE_NAME="lineageToPublicName.tsv"
# 
# echo "Downloading and decompressing consistent data files from ${RELEASE_DATE}..."
# 
# # Download and rename protobuf file
# wget -nc "${BASE_URL}/${PB_FILE_NAME}.gz" -O public-latest.all.masked.pb.gz
# gunzip -f public-latest.all.masked.pb.gz
# 
# # Download and rename metadata file
# wget -nc "${BASE_URL}/${METADATA_FILE_NAME}.gz" -O public-latest.metadata.tsv.gz
# gunzip -f public-latest.metadata.tsv.gz
# 
# # The reference genome is stable and does not need to be from a dated release
# wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/chromosomes/NC_045512v2.fa.gz
# gunzip -f NC_045512v2.fa.gz
# 
# # This is for information only since they are different than clades.tsv and samples.txt
# # wget -nc ${BASE_URL}/${CLADE_FILE_NAME}.gz
# # wget -nc ${BASE_URL}/${LINEAGE_FILE_NAME}.gz
# 
# matUtils summary --input-mat public-latest.all.masked.pb --clades clades.tsv
# matUtils summary --input-mat public-latest.all.masked.pb --samples samples.txt
# 
# echo "Data collection complete."

# @title
# Use menu item to fold this Cell (hide it from view)

# ==============================================================================
# Diagnostic Cell for Exploratory Analysis of Clade Sizes
# ==============================================================================
print("--- Analyzing Clade Size Distribution ---")

import os
import pandas as pd
import sys

# Reference file paths
SAMPLES_TSV = "samples.txt"
METADATA_FILE = "public-latest.metadata.tsv"

# Step 1: Load and parse sample list from the .pb file
if not os.path.isfile(SAMPLES_TSV):
    print(f"FATAL: Master sample list not found at: {SAMPLES_TSV}", file=sys.stderr)
    sys.exit(1)
master_sample_set = {line.strip().split('\t')[0] for line in open(SAMPLES_TSV) if line.strip()}

# Step 2: Load metadata and find common samples
metadata_df = pd.read_csv(
    METADATA_FILE, sep='\t', usecols=['strain', 'pangolin_lineage'], dtype=str
).dropna(subset=['strain', 'pangolin_lineage'])
filtered_metadata_df = metadata_df[metadata_df['strain'].isin(master_sample_set)]

if filtered_metadata_df.empty:
    print("FATAL: No common samples found.", file=sys.stderr)
    sys.exit(1)

# Step 3: Group by clade, count the members, and display the largest ones
clade_sizes = filtered_metadata_df.groupby('pangolin_lineage').size().sort_values(ascending=False)

print("\n--- Top 250 Largest Clades by Sample Count ---")
print(clade_sizes.head(250).to_string())
print("\n--- End of Analysis ---")

# This is based on top hits of the database only. Instead, use my master control file
# the top hits is fine to retrieve the master list, but then chunk the files and run Cell 4b instead

# This is the chunked file version, but run 4a or 4b Cell not both. 4a is to retrieve the control file only.

# ==============================================================================
# Cell 4: Subsampling Version (Robustly Handles Small Clades)
# ==============================================================================
print("--- Preparing a SUBSAMPLED set of batch files (Robust Version) ---")

import os
import pandas as pd
import sys

# --- Subsampling Parameters (You can adjust these) ---
RANGE_START = 201  # The starting rank (e.g., 1 is the largest)
RANGE_END = 250    # The ending rank (inclusive)
SAMPLES_PER_CLADE = 1000
# ----------------------------------------------------

# Reference file paths
SAMPLES_TSV = "samples.txt"
METADATA_FILE = "public-latest.metadata.tsv"
CLADES_CONTROL_FILE = "clades_to_process.txt"

print(f"Configuration: Selecting clades ranked {RANGE_START} through {RANGE_END} by size, and taking up to {SAMPLES_PER_CLADE} random samples from each.")

# Step 1: Load and parse sample list from the .pb file
master_sample_set = {line.strip().split('\t')[0] for line in open(SAMPLES_TSV) if line.strip()}
print(f"Loaded and parsed {len(master_sample_set)} unique strain IDs.")

# Step 2: Load metadata and find common samples
metadata_df = pd.read_csv(
    METADATA_FILE, sep='\t', usecols=['strain', 'pangolin_lineage'], dtype=str
).dropna(subset=['strain', 'pangolin_lineage'])
filtered_metadata_df = metadata_df[metadata_df['strain'].isin(master_sample_set)]
print(f"Found {len(filtered_metadata_df)} common samples between the files.")

if filtered_metadata_df.empty:
    print("FATAL: No common samples found.", file=sys.stderr)
    sys.exit(1)

# Step 3: Identify the clades within the specified rank range
grouped_metadata = filtered_metadata_df.groupby('pangolin_lineage')
# Sort all clades by size and select the specified range.
all_clade_sizes = grouped_metadata.size().sort_values(ascending=False)
# Slice using 0-based indexing. [RANGE_START - 1] gets the Nth item.
# Sort all clades by size to establish ranks
all_clade_sizes = grouped_metadata.size().sort_values(ascending=False)
# Slice the sorted list to get the desired range. Indices are 0-based, so RANGE_START must be adjusted.
major_clades = all_clade_sizes[RANGE_START - 1 : RANGE_END]
print(f"\nIdentified the following clades to sample from:\n{major_clades.to_string()}")

# Step 4: Generate subsampled per-clade lists and the control file
print("\nGenerating SUBSAMPLED per-clade sample lists...")
with open(CLADES_CONTROL_FILE, 'w') as f_control:
    for clade_name, group_df in grouped_metadata:
        if clade_name in major_clades.index:

            # Determine the number of samples to take
            num_samples_to_take = min(SAMPLES_PER_CLADE, len(group_df))

            # Sample using this safe, calculated number
            subsampled_group = group_df.sample(n=num_samples_to_take, random_state=42)

            safe_clade_name = "".join(c if c.isalnum() else "_" for c in clade_name)
            clade_sample_list_path = f"samples_{safe_clade_name}.txt"

            subsampled_group['strain'].to_csv(clade_sample_list_path, header=False, index=False)
            f_control.write(f"{safe_clade_name}\n")

print("\n--- Preparation of subsampled batch files is complete. ---")

# Run this to enable the persistent file option below.

# This and the Colab download file option depend on a connected runtime.
import os
from google.colab import drive

# Mount Google Drive
drive.mount('/content/drive')

# Confirm that the working directory has not changed
print(f"Current working directory remains: {os.getcwd()}")
print("Google Drive is mounted and accessible at /content/drive.")

# The following would be uncommented in Colab, but in its own Cell
#  of shell code only (not mixed with Python code).
# Commented out IPython magic to ensure Python compatibility.
# %%shell
# # ==============================================================================
# # Cell 5: Execute Batch Processing Loop (with Progress Indicator)
# # ==============================================================================
# set -e
# 
# # Define file paths
# PB_FILE="public-latest.all.masked.pb"
# REF_FASTA="NC_045512v2.fa"
# CLADES_CONTROL_FILE="clades_to_process.txt"
# FINAL_ALIGNED_FASTA="all_clades_aligned.fas"
# 
# # Ensure the final output file is empty before starting
# echo "Initializing final output file: $FINAL_ALIGNED_FASTA"
# # rm -f "$FINAL_ALIGNED_FASTA" # commented out in case of interruption of the Colab runtime
# touch "$FINAL_ALIGNED_FASTA"
# 
# if [ ! -f "$CLADES_CONTROL_FILE" ]; then
#     echo "FATAL: Clade control file not found at $CLADES_CONTROL_FILE"
#     exit 1
# fi
# 
# # --- MODIFICATION: Add progress tracking variables ---
# TOTAL_CLADES=$(wc -l < "$CLADES_CONTROL_FILE")
# clade_counter=0
# echo "Starting batch processing loop for ${TOTAL_CLADES} clades..."
# # ----------------------------------------------------
# 
# # Outer loop: Read each sanitized clade name
# while read -r safe_clade_name; do
#   # Supress the per-clade "Processing..." message to reduce verbosity
#   # echo "--- Processing Clade: ${safe_clade_name} ---"
# 
#   CLADE_SAMPLE_FILE="samples_${safe_clade_name}.txt"
#   TEMP_VCF="temp_${safe_clade_name}.vcf"
#   TEMP_VCF_GZ="${TEMP_VCF}.gz"
# 
#   if [ ! -s "$CLADE_SAMPLE_FILE" ]; then
#     # echo "Skipping ${safe_clade_name} due to empty or missing sample file."
#     continue
#   fi
# 
#   # Steps 1 & 2: Extract, compress, and index (output suppressed for cleaner log)
#   matUtils extract -i "$PB_FILE" --samples "$CLADE_SAMPLE_FILE" -v "$TEMP_VCF" >/dev/null 2>&1
#   bgzip --force "$TEMP_VCF"
#   bcftools index --force --tbi "$TEMP_VCF_GZ"
# 
#   # Step 3: Inner loop to generate consensus for each sample
#   while read -r sample_id; do
#     bcftools consensus -f "$REF_FASTA" --sample "$sample_id" "$TEMP_VCF_GZ" 2>/dev/null \
#       | sed "1s#.*#>${sample_id}|${safe_clade_name}#" >> "$FINAL_ALIGNED_FASTA"
#   done < "$CLADE_SAMPLE_FILE"
# 
#   # Step 4: Clean up intermediate files for this clade
#   rm "$TEMP_VCF_GZ" "$TEMP_VCF_GZ.tbi"
# 
#   # --- Print progress line ---
#   clade_counter=$((clade_counter + 1))
#   echo "[Progress ${clade_counter}/${TOTAL_CLADES}] Clade complete: ${safe_clade_name}"
#   # -----------------------------------------
# 
#   CHECKPOINT_INTERVAL=10
#   if (( clade_counter % CHECKPOINT_INTERVAL == 0 )); then
#     echo "--- Creating checkpoint at clade ${clade_counter} ---"
#     gzip -c "$FINAL_ALIGNED_FASTA" > "checkpoint_clade_${clade_counter}.fas.gz"
#     echo "--- Checkpoint file created: checkpoint_clade_${clade_counter}.fas.gz ---"
#   fi
# 
# done < "$CLADES_CONTROL_FILE"
# 
# # --- Post-Processing: Save to Google Drive Using Absolute Paths ---
# echo "--- Batch processing loop finished. ---"
# echo "Compressing and moving the final output file to Google Drive..."
# 
# # Define explicit, absolute paths to avoid ambiguity
# LOCAL_FILE_PATH="/content/${FINAL_ALIGNED_FASTA}"
# COMPRESSED_FILE_NAME="${FINAL_ALIGNED_FASTA}.gz"
# DRIVE_DESTINATION_PATH="/content/drive/MyDrive/"
# 
# # 1. Compress the final FASTA file
# gzip "${LOCAL_FILE_PATH}"
# 
# # 2. Move the resulting compressed file to Google Drive
# mv "/content/${COMPRESSED_FILE_NAME}" "${DRIVE_DESTINATION_PATH}"
# 
# echo "File successfully moved to ${DRIVE_DESTINATION_PATH}"
# echo "--- Workflow Complete ---"

# @title
# ==============================================================================
# (Standalone) Cell: Filter FASTA File to Remove an "Unassigned" Clade
# ==============================================================================
import os

# If running standalone without above steps, then mount google drive first as shown above

# --- Configuration ---
# Specify the input file generated by the previous steps.
INPUT_FASTA_PATH = "all_clades_aligned.fas"
# Specify the name for the new, filtered output file.
OUTPUT_FASTA_PATH = "all_clades_aligned_filtered.fas"
# -------------------

# Use lines below if the file is no longer available locally, but on google drive
!cp "/content/drive/MyDrive/{INPUT_FASTA_PATH}.gz" /content
!gunzip -f {INPUT_FASTA_PATH}

print(f"Starting filter process for: {INPUT_FASTA_PATH}")

# Verify that the input file exists before proceeding.
if not os.path.isfile(INPUT_FASTA_PATH):
    print(f"Error: Input file not found at '{INPUT_FASTA_PATH}'. Please ensure the file exists.")
else:
    sequences_read = 0
    sequences_written = 0

    with open(INPUT_FASTA_PATH, 'r') as f_in, open(OUTPUT_FASTA_PATH, 'w') as f_out:
        # A flag to determine if the current sequence block should be written.
        write_sequence = False

        for line in f_in:
            # Check if the line is a header line.
            if line.startswith('>'):
                sequences_read += 1
                # The clade is the last field when split by '|'.
                header_parts = line.strip().split('|')

                # Check if the header is well-formed and if the last part is 'Unassigned'.
                if len(header_parts) > 1 and header_parts[-1] == 'Unassigned':
                    # If it is an `Unassigned` clade, set the flag to False.
                    write_sequence = False
                else:
                    # Otherwise, set the flag to True and write this header.
                    write_sequence = True
                    f_out.write(line)
                    sequences_written += 1
            elif write_sequence:
                # If this is a sequence line and the flag is True, write the line.
                f_out.write(line)

    print("\n--- Filtering Complete ---")
    print(f"Total sequences read: {sequences_read}")
    print(f"Sequences removed (clade Unassigned): {sequences_read - sequences_written}")
    print(f"Sequences written to new file: {sequences_written}")
    print(f"Filtered data is available at: {OUTPUT_FASTA_PATH}")

# The following would be uncommented in Colab, but in its own Cell
#  of shell code only (not mixed with Python code).
# To copy the filtered output
# !gzip -f all_clades_aligned_filtered.fas
# !cp "/content/all_clades_aligned_filtered.fas.gz" /content/drive/MyDrive/

# @title
## Move Fasta File to Local Drive

# This and the Colab download file option depend on a connected runtime.

import os
from google.colab import drive

os.chdir('/content')
!gzip all_clades_aligned.fas

# Mount Google Drive
drive.mount('/content/drive')
os.chdir('/content/drive/MyDrive')
print(f"Current working directory: {os.getcwd()}")

# !mv /content/all_clades_aligned.fas.gz /content/drive/MyDrive/

# The following would be uncommented in Colab, but in its own Cell
#  of shell code only (not mixed with Python code).
# Commented out IPython magic to ensure Python compatibility.
# # @title
# %%shell
# 
# # This code not needed if the runtime is restarted.
# 
# # ==============================================================================
# # Cell 6: Verification and Final Cleanup
# # ==============================================================================
# set -e
# 
# # Define file paths using correct shell syntax (no spaces around '=')
# FINAL_ALIGNED_FASTA="all_clades_aligned.fas"
# CLADES_CONTROL_FILE="clades_to_process.txt"
# 
# echo "Verifying the final output file..."
# if [ ! -f "$FINAL_ALIGNED_FASTA" ]; then
#     echo "Error: Final output file was not created."
#     exit 1
# fi
# 
# # Count the number of sequences in the final FASTA file
# SEQUENCE_COUNT=$(grep -c ">" "$FINAL_ALIGNED_FASTA")
# echo "Total sequences in final file: $SEQUENCE_COUNT"
# 
# # Clean up the control file and all per-clade sample lists
# echo "Cleaning up control and sample list files..."
# rm "$CLADES_CONTROL_FILE"
# rm -f samples_*.txt
# 
# echo "--- Workflow Complete ---"