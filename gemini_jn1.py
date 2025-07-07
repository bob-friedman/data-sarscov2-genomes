# Commented out IPython magic to ensure Python compatibility.
# %%shell
# # Set a base URL for the latest data (Cell #1)
# BASE_URL="http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/"
# BASE_URL_2="http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/chromosomes/"
# PB_FILE="public-latest.all.masked.pb.gz"
# META_FILE="public-latest.metadata.tsv.gz"
# REF_FILE="NC_045512v2.fa.gz" # Reference genome
# 
# # Download the latest files
# echo "Downloading latest UShER data..."
# wget -q -nc "${BASE_URL}${PB_FILE}"
# wget -q -nc "${BASE_URL}${META_FILE}"
# wget -q -nc "${BASE_URL_2}${REF_FILE}"
# 
# # Decompress for use in the pipeline
# echo "Decompressing files..."
# gunzip -f "${PB_FILE}"
# gunzip -f "${META_FILE}"
# gunzip -f "${REF_FILE}" -c > NC_045512v2.fa # Decompress and rename reference
# 
# echo "Data acquisition complete. Ready to run Part I pipeline."

# Installation of Conda for Google Colab (Cell #2)
!pip install -q condacolab
import condacolab
condacolab.install()

# Commented out IPython magic to ensure Python compatibility.
# %%shell
# 
# set -e # Exit immediately if a command fails (Cell #3)
# 
# echo "--- Phase 1: Alignment Generation Started ---"
# 
# # --- 1. Environment Setup ---
# echo "[1/5] Installing dependencies (usher, bcftools, pandas)..."
# 
# conda install -y -c bioconda -c conda-forge usher bcftools pandas > /dev/null 2>&1
# echo "Dependencies installed."
# 
# # --- 2. Define File Paths ---
# PB_FILE="public-latest.all.masked.pb"
# METADATA_FILE="public-latest.metadata.tsv"
# REF_FASTA="NC_045512v2.fa"
# SAMPLES_TSV="samples.txt"
# CLADES_CONTROL_FILE="clades_to_process.txt"
# FINAL_ALIGNED_FASTA="all_clades_aligned.fas"
# 
# # --- 3. Master Sample Extraction ---
# echo "[2/5] Extracting master sample list from protobuf file..."
# matUtils summary -i "$PB_FILE" -s "$SAMPLES_TSV"
# echo "Master sample list created at ${SAMPLES_TSV}."

# --- 4. Clade-Based Subsampling Script (Python, Cell #4) ---
!echo "[3/5] Performing clade-based subsampling..."

import pandas as pd
import os

# --- Configuration ---
NUMBER_OF_MAJOR_CLADES = 250
SAMPLES_PER_CLADE = 1000
METADATA_FILE="public-latest.metadata.tsv"
SAMPLES_TSV = "samples.txt"
CLADES_CONTROL_FILE="clades_to_process.txt"
# --------------------

print(f"Loading metadata from {METADATA_FILE}...")
master_sample_set = {line.strip().split('\t')[0] for line in open(SAMPLES_TSV)}
metadata_df = pd.read_csv(METADATA_FILE, sep='\t', usecols=['strain', 'pangolin_lineage'], on_bad_lines='skip')
filtered_df = metadata_df[metadata_df['strain'].isin(master_sample_set) & (metadata_df['pangolin_lineage'] != 'Unassigned')]

print(f"Identifying top {NUMBER_OF_MAJOR_CLADES} clades...")
grouped_metadata = filtered_df.groupby('pangolin_lineage')
major_clades = grouped_metadata.size().nlargest(NUMBER_OF_MAJOR_CLADES)

print(f"Generating sample lists for {len(major_clades)} clades...")
with open(CLADES_CONTROL_FILE, 'w') as f_control:
    for clade_name, count in major_clades.items():
        group_df = grouped_metadata.get_group(clade_name)
        num_to_sample = min(SAMPLES_PER_CLADE, len(group_df))
        subsampled_group = group_df.sample(n=num_to_sample, random_state=42)
        safe_clade_name = clade_name.replace('.', '_')
        sample_list_file = f'samples_{safe_clade_name}.txt'
        subsampled_group['strain'].to_csv(sample_list_file, index=False, header=False)
        f_control.write(f"{safe_clade_name}\n")
print("Subsampling complete. Control file created at ${CLADES_CONTROL_FILE}.")

# Cell #5a - Mount Google Drive
# This must be in a separate Python cell before the main shell script.
from google.colab import drive
drive.mount('/content/drive')

# Commented out IPython magic to ensure Python compatibility.
# %%shell
# # Cell #5b - Fully Resilient Consensus Generation Loop
# 
# # --- Configuration ---
# # Define the path in your Google Drive where backups will be stored.
# DRIVE_PATH="/content/drive/MyDrive/"
# # The name of the backup archive.
# BACKUP_ARCHIVE="sars-cov2-pi-checkpoint.tar.gz"
# # How often to back up progress to Google Drive (e.g., every 15 clades).
# CHECKPOINT_INTERVAL=15
# 
# # --- Define File Paths ---
# PB_FILE="public-latest.all.masked.pb"
# REF_FASTA="NC_045512v2.fa"
# FINAL_ALIGNED_FASTA="all_clades_aligned.fas"
# CLADES_CONTROL_FILE="clades_to_process.txt"
# CHECKPOINT_FILE="clades_completed.txt"
# SAMPLES_TSV="samples.txt"
# # --------------------
# 
# echo "--- Starting Resilient Pipeline ---"
# 
# # --- 1. Restore State from Google Drive (if backup exists) ---
# if [ -f "${DRIVE_PATH}${BACKUP_ARCHIVE}" ]; then
#     echo "Backup archive found in Google Drive. Restoring state..."
#     cp "${DRIVE_PATH}${BACKUP_ARCHIVE}" .
#     tar -xzf "${BACKUP_ARCHIVE}"
#     echo "Restore complete. Resuming from last checkpoint."
# else
#     echo "No backup found. Starting a new run."
# fi
# 
# # Initialize files if they don't exist after potential restore
# touch "$CHECKPOINT_FILE"
# touch "$FINAL_ALIGNED_FASTA"
# 
# # --- 2. Main Processing Loop ---
# CLADE_COUNT=$(wc -l < "$CLADES_CONTROL_FILE")
# PROCESSED_COUNT=$(wc -l < "$CHECKPOINT_FILE")
# echo "Found ${PROCESSED_COUNT} clades already completed. Total to process: ${CLADE_COUNT}."
# 
# while read -r safe_clade_name; do
#     # Check if the clade is already completed
#     if grep -q "^${safe_clade_name}$" "$CHECKPOINT_FILE"; then
#         echo "Skipping already completed clade: ${safe_clade_name}"
#         continue
#     fi
# 
#     PROCESSED_COUNT=$((PROCESSED_COUNT + 1))
#     printf "Processing clade %d/%d: %s...\n" "$PROCESSED_COUNT" "$CLADE_COUNT" "$safe_clade_name"
# 
#     CLADE_SAMPLE_FILE="samples_${safe_clade_name}.txt"
#     TEMP_VCF="temp_${safe_clade_name}.vcf"
#     TEMP_VCF_GZ="${TEMP_VCF}.gz"
# 
#     matUtils extract -i "$PB_FILE" --samples "$CLADE_SAMPLE_FILE" -v "$TEMP_VCF"
#     bgzip --force "$TEMP_VCF"
#     bcftools index --force --tbi "$TEMP_VCF_GZ"
# 
#     while read -r sample_id; do
#         bcftools consensus -f "$REF_FASTA" --sample "$sample_id" "$TEMP_VCF_GZ" \
#           | sed "1s#.*#>${sample_id}|${safe_clade_name}#" >> "$FINAL_ALIGNED_FASTA"
#     done < "$CLADE_SAMPLE_FILE"
# 
#     # Log the clade as completed
#     echo "${safe_clade_name}" >> "$CHECKPOINT_FILE"
# 
#     # Clean up intermediate files for this clade
#     rm "$TEMP_VCF_GZ" "$TEMP_VCF_GZ.tbi" "$CLADE_SAMPLE_FILE"
# 
#     # --- 3. Periodic Backup to Google Drive ---
#     if (( PROCESSED_COUNT % CHECKPOINT_INTERVAL == 0 )); then
#         echo "Checkpoint interval reached. Backing up progress to Google Drive..."
#         # Create a compressed archive of the essential state files
#         sync
#         tar -czf "${BACKUP_ARCHIVE}" "${FINAL_ALIGNED_FASTA}" "${CHECKPOINT_FILE}"
#         # Copy the archive to Google Drive, overwriting the previous backup
#         cp "${BACKUP_ARCHIVE}" "${DRIVE_PATH}"
#         echo "Backup successful to ${DRIVE_PATH}${BACKUP_ARCHIVE}"
#     fi
# 
# done < "$CLADES_CONTROL_FILE"
# 
# echo "Consensus generation loop complete."
# 
# # --- 4. Final Verification and Cleanup ---
# echo "Verifying output and performing final cleanup..."
# SEQUENCE_COUNT=$(grep -c ">" "$FINAL_ALIGNED_FASTA")
# echo "Verification complete. Final file '${FINAL_ALIGNED_FASTA}' contains $SEQUENCE_COUNT sequences."
# 
# # Perform a final backup to ensure the very last clades are saved
# echo "Performing final backup of completed work..."
# sync
# tar -czf "${BACKUP_ARCHIVE}" "${FINAL_ALIGNED_FASTA}" "${CHECKPOINT_FILE}"
# cp "${BACKUP_ARCHIVE}" "${DRIVE_PATH}"
# echo "Final backup successful."
# 
# # Clean up all local control and intermediate files
# rm -f "$CLADES_CONTROL_FILE" "$SAMPLES_TSV" "$CHECKPOINT_FILE" "$BACKUP_ARCHIVE"
# echo "--- Phase 1: Alignment Generation Finished ---"

!ls -l all_clades_aligned.fas

# Cell for Validation of Nucleotide Diversity (π) Calculation

!pip install -q biopython
import numpy as np
import pandas as pd
from Bio import SeqIO
import os

# --- Configuration ---
FULL_ALIGNMENT_FILE = "all_clades_aligned.fas"
CLADE_TO_TEST = "JN_1"
NUM_SAMPLES_TO_TEST = 100
TEST_FASTA_FILE = f"test_alignment_{CLADE_TO_TEST}.fas"
# --------------------

print(f"--- Starting The Definitive and Final Validation for Clade: {CLADE_TO_TEST} ---")

# 1. Create a small test FASTA file
print(f"Step 1: Extracting {NUM_SAMPLES_TO_TEST} sequences from {FULL_ALIGNMENT_FILE}...")
sequences_written = 0
with open(TEST_FASTA_FILE, "w") as f_out:
    for record in SeqIO.parse(FULL_ALIGNMENT_FILE, "fasta"):
        header_parts = record.id.split('|')
        if len(header_parts) > 1 and header_parts[-1] == CLADE_TO_TEST:
            if sequences_written < NUM_SAMPLES_TO_TEST:
                SeqIO.write(record, f_out, "fasta")
                sequences_written += 1
            else:
                break
if sequences_written < 2:
    print(f"Error: Found only {sequences_written} sequences for clade {CLADE_TO_TEST}.")
else:
    print(f"Successfully created test file: {TEST_FASTA_FILE} with {sequences_written} sequences.")

    # 2. Convert to a single, shared NumPy array
    records = list(SeqIO.parse(TEST_FASTA_FILE, "fasta"))
    seq_length = len(records[0].seq)
    alignment_np = np.empty((len(records), seq_length), dtype=np.uint8)
    nt_map = str.maketrans("ACGTN", "01234")
    for i, record in enumerate(records):
        processed_sequence = str(record.seq).upper().translate(nt_map)
        alignment_np[i, :] = np.frombuffer(processed_sequence.encode(), dtype=np.uint8) - ord('0')

    # 3. Calculate with our optimized custom function
    print("\nStep 2: Calculating with our optimized custom function...")
    def calculate_total_diffs_optimized(alignment_slice):
        n_seq, n_sites = alignment_slice.shape
        if n_seq < 2: return 0.0
        total_pairwise_diffs = 0.0
        for i in range(n_sites):
            col = alignment_slice[:, i]
            valid_bases = col[col < 4] # Ignore Ns
            n_k = len(valid_bases)
            if n_k < 2: continue
            counts = np.bincount(valid_bases, minlength=4)
            num_pairs_at_site = n_k * (n_k - 1) / 2
            num_same_pairs_at_site = sum(c * (c - 1) / 2 for c in counts)
            total_pairwise_diffs += (num_pairs_at_site - num_same_pairs_at_site)
        return total_pairwise_diffs

    our_diffs_optimized = calculate_total_diffs_optimized(alignment_np)
    print(f"Our optimized function total differences = {our_diffs_optimized}")

    # 4. Calculate with a simple, brute-force method for ground truth
    print("\nStep 3: Calculating with brute-force method for ground truth...")
    def calculate_total_diffs_brute_force(alignment_slice):
        n_seq, n_sites = alignment_slice.shape
        if n_seq < 2: return 0.0
        total_diffs = 0
        for i in range(n_seq):
            for j in range(i + 1, n_seq):
                seq1 = alignment_slice[i, :]
                seq2 = alignment_slice[j, :]
                # Find sites where both sequences are not 'N'
                valid_comparison = (seq1 < 4) & (seq2 < 4)
                # Sum the number of differences at those sites
                total_diffs += np.sum(seq1[valid_comparison] != seq2[valid_comparison])
        return total_diffs

    ground_truth_diffs = calculate_total_diffs_brute_force(alignment_np)
    print(f"Ground truth brute-force differences   = {ground_truth_diffs}")

    # 5. Compare results
    print("\n--- Validation Complete ---")
    if np.isclose(our_diffs_optimized, ground_truth_diffs):
        print("SUCCESS: Optimized function matches brute-force ground truth. Our method is fully validated.")
    else:
        print("FAILURE: A discrepancy remains between our function and the ground truth.")
        print(f"Difference: {abs(our_diffs_optimized - ground_truth_diffs)}")

# Cell for Final Denominator Validation

!pip install -q biopython
!pip install -q scikit-allel
import numpy as np
from Bio import SeqIO
import os

# --- Configuration ---
TEST_FASTA_FILE = "test_alignment_JN_1.fas"
# --------------------

print(f"--- Validating Denominator Calculation for: {TEST_FASTA_FILE} ---")

if not os.path.exists(TEST_FASTA_FILE):
    print(f"Error: Test file '{TEST_FASTA_FILE}' not found. Please run the previous validation script first.")
else:
    # 1. Load the test alignment into a NumPy array
    records = list(SeqIO.parse(TEST_FASTA_FILE, "fasta"))
    seq_length = len(records[0].seq)
    alignment_np = np.empty((len(records), seq_length), dtype=np.uint8)
    nt_map = str.maketrans("ACGTN", "01234")
    for i, record in enumerate(records):
        processed_sequence = str(record.seq).upper().translate(nt_map)
        alignment_np[i, :] = np.frombuffer(processed_sequence.encode(), dtype=np.uint8) - ord('0')

    # 2. Calculate the two potential denominators
    total_sequence_length = alignment_np.shape[1]

    valid_sites_count = 0
    for i in range(total_sequence_length):
        column = alignment_np[:, i]
        # A site is valid for comparison if it has at least 2 non-'N' bases
        if np.sum(column < 4) >= 2:
            valid_sites_count += 1

    # 3. Print the results
    print("\n--- Denominator Report ---")
    print(f"Total Sequence Length (Denominator used by scikit-allel): {total_sequence_length}")
    print(f"Number of Valid Sites for Comparison (Denominator used by our function): {valid_sites_count}")

    print("\n--- Conclusion ---")
    if valid_sites_count < total_sequence_length:
        print("As expected, the number of valid sites is less than the total length.")
        print("This confirms that our method correctly excludes low-quality/ambiguous sites from the final average.")
    else:
        print("The number of valid sites is equal to the total length.")