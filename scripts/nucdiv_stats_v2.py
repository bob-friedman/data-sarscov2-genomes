# @title
'''
# Part 1: Revised Standalone Verification and Filtering Script

This script has been modified to produce a third artifact: `verified_sequence_order.txt`.
This file contains the sequence IDs in the exact order they appear in the clean FASTA file, serving as
an essential index for aligning metadata with the sequence data in the next stage.
'''

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
# #### **Standalone Verification and Filtering Cell**
# 
# # @title Standalone FASTA Verification and Filtering
# %%shell
# # This block handles file acquisition and assembly on the Colab local disk.
# 
# # --- Configuration ---
# DRIVE_PATH="/content/drive/MyDrive/"
# # Add "part" to filename if absent
# FILE_PARTS=(
#   "all_clades_aligned_part1.fas.gz"
#   "all_clades_aligned_part2.fas.gz"
#   "all_clades_aligned_part3.fas.gz"
#   "all_clades_aligned_part4.fas.gz"
#   "all_clades_aligned_part5.fas.gz"
# )
# RAW_FASTA_NAME="all_clades_aligned_combined.fas"
# # --------------------
# 
# rm -f "${RAW_FASTA_NAME}"
# echo "--- Starting File Assembly on Colab Local Disk ---"
# 
# for part in "${FILE_PARTS[@]}"; do
#   if [ -f "${DRIVE_PATH}${part}" ]; then
#     echo "Copying ${part} from Google Drive..."
#     cp "${DRIVE_PATH}${part}" .
#   else
#     echo "Warning: File not found in Google Drive: ${DRIVE_PATH}${part}"
#   fi
# done
# 
# echo "Decompressing and concatenating parts..."
# gunzip -c *.fas.gz > "${RAW_FASTA_NAME}"
# rm -f *.fas.gz # Clean up intermediate parts
# echo "--- File Assembly Complete. Raw file is at /content/${RAW_FASTA_NAME} ---"

# This Python block handles the memory-efficient filtering and index creation.
import os
import csv

# --- Configuration ---
INPUT_FASTA = "all_clades_aligned_combined.fas"
CLEAN_FASTA_OUTPUT = "verified_sequences.fas"
LOG_FILE_OUTPUT = "verification_log.csv"
ORDER_INDEX_OUTPUT = "verified_sequence_order.txt" # New output file
# --------------------

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
    """
    Streams a FASTA file, verifies each sequence, writes valid ones to a new
    file, logs discards, and creates an ordered index of valid IDs.
    """
    print(f"Starting verification of '{input_path}'...")
    if not os.path.isfile(input_path):
        print(f"FATAL: Input file not found.")
        return

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
                if ref_len is None:
                    ref_len = len(sequence)
                elif len(sequence) != ref_len:
                    error_reason = f"Incorrect length (is {len(sequence)}, expected {ref_len})"

            if error_reason:
                errors += 1
                log_writer.writerow([seq_id, error_reason])
            else:
                written += 1
                f_clean.write(f">{header}\n{sequence}\n")
                f_order.write(f"{seq_id}\n")

    print("\n--- Verification and Filtering Complete ---")
    print(f"Total sequences processed: {processed}")
    print(f"Sequences written to new file: {written}")
    print(f"Sequences discarded (errors): {errors}")
    print(f"Clean data: {clean_path}")
    print(f"Error log: {log_path}")
    print(f"Ordered ID index: {order_path}")

# --- Main Execution ---
verify_and_filter_fasta(INPUT_FASTA, CLEAN_FASTA_OUTPUT, LOG_FILE_OUTPUT, ORDER_INDEX_OUTPUT)

# The following would be uncommented in Colab, but in its own Cell
#  of shell code only (not mixed with Python code).
# Commented out IPython magic to ensure Python compatibility.
# %%shell
# # This block compresses the final outputs and moves all three to Google Drive.
# 
# CLEAN_FASTA="verified_sequences.fas"
# LOG_FILE="verification_log.csv"
# ORDER_FILE="verified_sequence_order.txt"
# DRIVE_PATH="/content/drive/MyDrive/"
# 
# echo "--- Compressing and Moving Final Artifacts to Google Drive ---"
# # compression finishes in ~15 minutes
# gzip -f "${CLEAN_FASTA}"
# mv "${CLEAN_FASTA}.gz" "${DRIVE_PATH}"
# mv "${LOG_FILE}" "${DRIVE_PATH}"
# mv "${ORDER_FILE}" "${DRIVE_PATH}"
# echo "Successfully moved artifacts to ${DRIVE_PATH}"
# rm -f "all_clades_aligned_combined.fas"

# @title
'''
# Part 2: Revised Partitioned Nucleotide Diversity Script

This script is now designed for partitioned analysis. It uses the three files generated above to calculate
nucleotide diversity for user-defined groups (e.g., by clade, time, and region) in a memory-efficient manner.
'''

# The following would be uncommented in Colab, but in its own Cell
#  of shell code only (not mixed with Python code).
# Commented out IPython magic to ensure Python compatibility.
# # @title
# # (Optional) Download metadata file if needed
# 
# %%shell
# RELEASE_DATE="2025/06/29"
# BASE_URL="http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/${RELEASE_DATE}/"
# METADATA_FILE_NAME="public-${RELEASE_DATE//\//-}.metadata.tsv"
# 
# echo "Downloading and decompressing consistent data files from ${RELEASE_DATE}..."
# 
# # Download and rename metadata file
# wget -nc "${BASE_URL}/${METADATA_FILE_NAME}.gz" -O public-latest.metadata.tsv.gz
# gunzip -f public-latest.metadata.tsv.gz
# 
# echo "Data collection complete."

# Cell 1: Setup and Data Acquisition

# @title Partitioned Nucleotide Diversity (π) - Setup
import numpy as np
import pandas as pd
import os

# --- Configuration ---
DRIVE_PATH = "/content/drive/MyDrive/"
# Inputs from verification step
CLEAN_FASTA_GZ = "verified_sequences.fas.gz"
ORDER_INDEX_FILE = "verified_sequence_order.txt"
# This is your main, original metadata file
METADATA_TSV = "public-latest.metadata.tsv"

# Local file names
CLEAN_FASTA_LOCAL = "verified_sequences.fas"
NUMPY_BINARY_FILE = "alignment.mmap"
# --------------------

# Acquire all necessary files from Google Drive
print("Acquiring input files from Google Drive...")
!cp "{DRIVE_PATH}{CLEAN_FASTA_GZ}" .
!cp "{DRIVE_PATH}{ORDER_INDEX_FILE}" .
!cp "{DRIVE_PATH}{METADATA_TSV}" .
!gunzip -f "{CLEAN_FASTA_GZ}"
print("All data ready for processing.")

# Cell 2: Pre-processing (FASTA to Binary)

'''
This cell is functionally unchanged but is a required step. It converts the large FASTA alignment
into a NumPy binary file for efficient access.
'''

# @title Partitioned Nucleotide Diversity (π) - Pre-processing

def fasta_to_numpy_binary(fasta_path, output_path):
    """Converts a FASTA alignment to a memory-mapped NumPy binary file."""
    print("Starting conversion from FASTA to memory-mapped binary format...")
    num_sequences, seq_length = 0, 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'): num_sequences += 1
            elif not seq_length: seq_length = len(line.strip())

    if num_sequences == 0: return None, 0, 0
    print(f"Alignment dimensions: {num_sequences} sequences, {seq_length} sites.")

    mmap_array = np.memmap(output_path, dtype=np.uint8, mode='w+', shape=(num_sequences, seq_length))
    nt_map = str.maketrans("ACGTN", "01234")

    current_seq_index = 0
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith('>'): continue
            # previous line of code that works, but shows warning logs about deprecated method
            # mmap_array[current_seq_index, :] = np.fromstring(line.strip().upper().translate(nt_map), dtype=np.uint8)
            # updated with two lines of code below to replace above line
            processed_sequence = line.strip().upper().translate(nt_map)
            mmap_array[current_seq_index, :] = np.frombuffer(processed_sequence.encode(), dtype=np.uint8) - ord('0')
            current_seq_index += 1

    mmap_array.flush()
    print(f"Successfully created binary file at: {output_path}")
    return output_path, num_sequences, seq_length

binary_file, n_seq, n_len = fasta_to_numpy_binary(CLEAN_FASTA_LOCAL, NUMPY_BINARY_FILE)

import pandas as pd
import gc

# --- Configuration ---
ORDER_INDEX_FILE = "verified_sequence_order.txt"
METADATA_TSV = "public-latest.metadata.tsv"
CHUNK_SIZE = 500000 # Process 500,000 rows at a time. Adjust if needed.
# ---------------------

print("--- Aligning Metadata with Sequence Data (Memory-Efficient) ---")

# 1. Load sequence IDs into a set for fast lookups.
ordered_ids_set = set(pd.read_csv(ORDER_INDEX_FILE, header=None, names=['strain'])['strain'])

# 2. Process the large metadata file in chunks to avoid loading it all at once.
cols_to_use = ['strain', 'pangolin_lineage', 'date', 'country']
col_dtypes = {'strain': 'object', 'pangolin_lineage': 'category', 'date': 'object', 'country': 'category'}

filtered_chunks = []
print(f"Reading '{METADATA_TSV}' in chunks of {CHUNK_SIZE} rows...")

with pd.read_csv(METADATA_TSV, sep='\t', usecols=cols_to_use, dtype=col_dtypes, chunksize=CHUNK_SIZE) as reader:
    for chunk in reader:
        # Normalize the 'strain' column for the current chunk.
        chunk['strain'] = chunk['strain'].str.split('|').str[0]

        # Filter the chunk to keep only rows whose strain is in our set.
        filtered_chunk = chunk[chunk['strain'].isin(ordered_ids_set)]
        filtered_chunks.append(filtered_chunk)
        print(f"  - Processed a chunk, found {len(filtered_chunk)} matching records.")

# 3. Concatenate the filtered chunks into a single DataFrame.
print("Concatenating filtered chunks...")
df_meta_filtered = pd.concat(filtered_chunks, ignore_index=True)

# Clean up intermediate list of chunks.
del filtered_chunks
gc.collect()

# --- De-duplicate the metadata ---
# The metadata file has multiple entries for some strains. We keep only the first one.
print(f"De-duplicating metadata. Shape before: {df_meta_filtered.shape}")
df_meta_filtered.drop_duplicates(subset=['strain'], keep='first', inplace=True)
print(f"Shape after de-duplication: {df_meta_filtered.shape}")
# ------------------------------------------

# 4. Reorder the filtered metadata to match the sequence file order.
# Note: We now need the original ordered list, not the set.
ordered_ids_list = pd.read_csv(ORDER_INDEX_FILE, header=None, names=['strain'])['strain'].tolist()
df_meta_aligned = df_meta_filtered.set_index('strain').loc[ordered_ids_list].reset_index()

# Clean up to save memory.
del df_meta_filtered
del ordered_ids_list
gc.collect()

print(f"Metadata aligned. Shape of final metadata table: {df_meta_aligned.shape}")

# 5. Define analytical partitions
print("\n--- Defining Analytical Partitions ---")
df_meta_aligned['date'] = pd.to_datetime(df_meta_aligned['date'], errors='coerce')
df_meta_aligned.dropna(subset=['date'], inplace=True)
df_meta_aligned['time_bin'] = df_meta_aligned['date'].dt.to_period('2W')
grouped = df_meta_aligned.groupby(['pangolin_lineage', 'time_bin'])
print(f"Created {len(grouped)} groups for analysis.")

# Cell 4: Partitioned Calculation and Results

'''
This cell iterates through the defined groups, calculates π for each one, and compiles the results into a final table.
'''

import numpy as np
import pandas as pd

# --- Configuration ---
# Define the path to the binary file and its dimensions from the previous step.
binary_file = "alignment.mmap"
n_seq = 248103  # Number of sequences from your previous logs
n_len = 29903   # Number of sites (sequence length) from your previous logs
# ---------------------

# @title Partitioned Nucleotide Diversity (π) - Calculation
def calculate_pi_robust(alignment_slice):
    """
    Calculates nucleotide diversity (π) using a validated, robust, per-site method.
    This function correctly handles missing data ('N's).
    """
    n_seq, n_sites = alignment_slice.shape
    if n_seq < 2:
        return 0.0

    pi_sum_at_valid_sites = 0.0
    sites_compared = 0

    # Iterate over each site (column) in the alignment
    for i in range(n_sites):
        col = alignment_slice[:, i]
        # Consider only valid bases (A,C,G,T), ignoring 'N's
        valid_bases = col[col < 4]
        n_k = len(valid_bases)

        # A comparison is only possible if there are at least 2 valid sequences
        if n_k < 2:
            continue
        
        sites_compared += 1
        
        # Calculate allele frequencies at this site
        counts = np.bincount(valid_bases, minlength=4)
        
        # Calculate π for this site using Tajima's estimator: (n/(n-1)) * (1 - sum(p_i^2))
        sum_freq_sq = np.sum((counts / n_k)**2)
        pi_site = (n_k / (n_k - 1)) * (1 - sum_freq_sq)
        
        pi_sum_at_valid_sites += pi_site

    # Average the per-site π values over the number of valid sites
    final_pi = pi_sum_at_valid_sites / sites_compared if sites_compared > 0 else 0.0
    return final_pi

# --- Main Calculation Loop ---
print("\n--- Starting Partitioned π Calculation ---")
# Open the binary file once in read-only mode
alignment = np.memmap(binary_file, dtype=np.uint8, mode='r', shape=(n_seq, n_len))

results = []
group_counter = 0
for name, group_df in grouped:
    group_counter += 1
    if len(group_df) < 2: # Cannot calculate diversity on less than 2 samples
        continue

    # Get the integer indices for this group from the aligned metadata DataFrame
    indices = group_df.index.tolist()

    # Use these indices to select the correct rows from the memory-mapped array
    alignment_slice = alignment[indices, :]

    pi = calculate_pi_robust(alignment_slice)

    # Add this line to print the result for each group immediately.
    print(f"L: {name[0]}, T: {str(name[1])}, Pi: {pi:.6f}, N: {len(group_df)}")

    results.append({
        'lineage': name[0],
        'time_bin': str(name[1]), # Convert period to string for CSV
        'pi': pi,
        'n_samples': len(group_df)
    })

    if group_counter % 100 == 0:
        print(f"  Processed {group_counter}/{len(grouped)} groups...")

print("\n--- Calculation Complete ---")
results_df = pd.DataFrame(results)
print("Sample of results:")
print(results_df.head())

# --- Save Results and Cleanup ---
results_df.to_csv("nucleotide_diversity_results.csv", index=False)
print("\nResults saved to nucleotide_diversity_results.csv")
!mv nucleotide_diversity_results.csv "{DRIVE_PATH}"
# do not remove these files if a diagnostic below is needed
# !rm -f "{CLEAN_FASTA_LOCAL}" "{NUMPY_BINARY_FILE}" "{ORDER_INDEX_FILE}" "{METADATA_TSV}"
print("Cleanup of local files complete.")

# @title Standalone Diagnostic Cell for High π Values (Corrected)

import numpy as np
import pandas as pd

# --- Configuration ---
PI_THRESHOLD = 0.001
RESULTS_FILE = "nucleotide_diversity_results.csv"
BINARY_FILE = "alignment.mmap"
N_SEQ = 248103
N_LEN = 29903
# --------------------

print(f"--- Diagnosing Groups with π > {PI_THRESHOLD} ---")

# 1. Load the results and find the high-pi groups.
results_df = pd.read_csv(RESULTS_FILE)
high_pi_groups = results_df[results_df['pi'] > PI_THRESHOLD]

if high_pi_groups.empty:
    print("No groups found exceeding the pi threshold.")
else:
    # 2. Re-create the grouped object to get indices.
    # This logic must exactly match the main calculation cell.
    df_meta_aligned['date'] = pd.to_datetime(df_meta_aligned['date'], errors='coerce')
    df_meta_aligned.dropna(subset=['date'], inplace=True)
    df_meta_aligned['time_bin'] = df_meta_aligned['date'].dt.to_period('2W')
    grouped = df_meta_aligned.groupby(['pangolin_lineage', 'time_bin'])

    # 3. Open the alignment file.
    alignment = np.memmap(BINARY_FILE, dtype=np.uint8, mode='r', shape=(n_seq, n_len))

    # 4. Iterate through only the high-pi groups for diagnosis.
    for index, row in high_pi_groups.iterrows():
        lineage = row['lineage']

        # Correctly reconstruct the Period object with '2W' frequency.
        time_bin_period = pd.Period(row['time_bin'].split('/')[0], freq='2W')
        group_name = (lineage, time_bin_period)

        try:
            group_df = grouped.get_group(group_name)
            indices = group_df.index.tolist()
            alignment_slice = alignment[indices, :]

            max_diff = 0
            num_samples = alignment_slice.shape[0]
            outlier_pair = (None, None)

            for i in range(num_samples):
                for j in range(i + 1, num_samples):
                    valid_comparison = (alignment_slice[i, :] < 4) & (alignment_slice[j, :] < 4)
                    diff = np.sum(np.not_equal(alignment_slice[i, valid_comparison], alignment_slice[j, valid_comparison]))

                    if diff > max_diff:
                        max_diff = diff
                        outlier_pair = (group_df.iloc[i]['strain'], group_df.iloc[j]['strain'])

            print(f"\n[High Pi Detected] Lineage: {lineage}, Time Bin: {row['time_bin']}, Pi: {row['pi']:.6f}")
            print(f"  - Samples in group: {num_samples}")
            print(f"  - Maximum difference between any two sequences: {max_diff} sites")
            print(f"  - Most divergent pair: {outlier_pair[0]} and {outlier_pair[1]}")

        except KeyError:
            # This warning should now only appear for legitimate data inconsistencies.
            print(f"\nWarning: Could not find group {group_name} in the metadata. Skipping diagnosis.")

# @title Filter Final Results for Robustness

import pandas as pd

# --- Configuration ---
RESULTS_FILE = "nucleotide_diversity_results.csv"
FILTERED_RESULTS_FILE = "nucleotide_diversity_results_filtered.csv"
DRIVE_PATH = "/content/drive/MyDrive/"

# Define filtering thresholds.
# A minimum of 10 samples is a common threshold for statistical reliability.
MIN_SAMPLES_THRESHOLD = 10
# The threshold for high pi values, consistent with the diagnostic script.
PI_THRESHOLD = 0.001
# --------------------

# Load the results
results_df = pd.read_csv(RESULTS_FILE)
initial_rows = len(results_df)
print(f"Loaded initial results with {initial_rows} rows.")

# Apply the filters
filtered_df = results_df[
    (results_df['n_samples'] >= MIN_SAMPLES_THRESHOLD) &
    (results_df['pi'] <= PI_THRESHOLD)
].copy()

final_rows = len(filtered_df)
rows_removed = initial_rows - final_rows
print(f"Applied filters (n_samples >= {MIN_SAMPLES_THRESHOLD}, pi <= {PI_THRESHOLD}).")
print(f"Removed {rows_removed} rows. Final dataset contains {final_rows} rows.")

# Save the filtered results
filtered_df.to_csv(FILTERED_RESULTS_FILE, index=False)
print(f"\nFiltered results saved to: {FILTERED_RESULTS_FILE}")

# Move the final, clean file to Google Drive
!mv "{FILTERED_RESULTS_FILE}" "{DRIVE_PATH}"
print(f"Successfully moved filtered results to {DRIVE_PATH}")

# @title Generate Heatmap of Nucleotide Diversity by Lineage and Time

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
FILTERED_RESULTS_FILE = "nucleotide_diversity_results_filtered.csv"
# Select the top N lineages to display for clarity.
NUM_LINEAGES_TO_DISPLAY = 20
# --------------------

# 1. Load the filtered data.
pi_results = pd.read_csv(FILTERED_RESULTS_FILE)

# Convert time_bin strings to Period objects for correct sorting.
pi_results['time_bin'] = pd.to_datetime(pi_results['time_bin'].str.split('/').str[0]).dt.to_period('2W')

# 2. Identify the top N lineages based on total sample count.
top_lineages = pi_results.groupby('lineage')['n_samples'].sum().nlargest(NUM_LINEAGES_TO_DISPLAY).index
pi_results_subset = pi_results[pi_results['lineage'].isin(top_lineages)]

# 3. Create the pivot table for the heatmap.
pivot_df = pi_results_subset.pivot_table(
    index='lineage',
    columns='time_bin',
    values='pi',
    aggfunc='mean'
)
# Ensure columns (time bins) are sorted chronologically.
pivot_df = pivot_df.reindex(sorted(pivot_df.columns), axis=1)

# 4. Generate the heatmap.
plt.figure(figsize=(24, 14))
heatmap = sns.heatmap(
    pivot_df,
    annot=False,  # Annotation is impractical for this many cells.
    cmap='viridis',
    linewidths=.5,
    cbar_kws={'label': 'Nucleotide Diversity (π)'}
)
plt.title(f'Nucleotide Diversity (π) for the Top {NUM_LINEAGES_TO_DISPLAY} Lineages Over Time', fontsize=16)
plt.xlabel('Time Bin (Bi-weekly)', fontsize=12)
plt.ylabel('Pango Lineage', fontsize=12)

# Improve x-axis label readability.
# Convert Period objects back to strings for display.
x_labels = [f"{period.start_time.strftime('%Y-%m-%d')}" for period in pivot_df.columns]
# Display only every nth label to avoid overcrowding.
n = 4
plt.xticks(ticks=range(len(x_labels))[::n], labels=x_labels[::n], rotation=45, ha='right')

plt.tight_layout()
plt.show()
