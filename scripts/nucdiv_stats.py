# This cell will copy our newly generated alignment from Google Drive to the local Colab environment and decompress it.

# Cell 1: Setup and Data Acquisition
from google.colab import drive
import os

print("--- Phase 2: Setup and Data Acquisition ---")

# Mount Google Drive
drive.mount('/content/drive')

# --- Configuration ---
DRIVE_PATH = "/content/drive/MyDrive/"
# This is the final backup archive from our Part I pipeline
BACKUP_ARCHIVE = "sars-cov2-pi-checkpoint.tar.gz"

# Local file names
RAW_FASTA_FROM_ARCHIVE = "all_clades_aligned.fas"
METADATA_TSV_LOCAL = "public-latest.metadata.tsv" # Assumes this is in Colab environment
# --------------------

# Restore the final alignment from the backup archive
if os.path.exists(f"{DRIVE_PATH}{BACKUP_ARCHIVE}"):
    print(f"Restoring '{RAW_FASTA_FROM_ARCHIVE}' from backup archive in Google Drive...")
    !cp "{DRIVE_PATH}{BACKUP_ARCHIVE}" .
    !tar -xzf "{BACKUP_ARCHIVE}"
    # We only need the FASTA file from the archive for this part
    if os.path.exists(RAW_FASTA_FROM_ARCHIVE):
        print("Successfully restored the final alignment file.")
    else:
        print("FATAL: Archive was found but did not contain the alignment file.")
else:
    print("FATAL: Backup archive not found in Google Drive.")

# This cell runs the verification script on our 7.6 GB alignment. It will produce the clean FASTA file and the
# crucial ordered index file. This step may take 10-15 minutes.

# Cell 2: Standalone Verification and Filtering
!pip install -q biopython
import os
import csv
from Bio import SeqIO

print("\n--- Starting Verification of Final Alignment ---")

# --- Configuration ---
INPUT_FASTA = "all_clades_aligned.fas"
CLEAN_FASTA_OUTPUT = "verified_sequences.fas"
LOG_FILE_OUTPUT = "verification_log.csv"
ORDER_INDEX_OUTPUT = "verified_sequence_order.txt"
# --------------------

def verify_and_filter_fasta(input_path, clean_path, log_path, order_path):
    processed, written, errors = 0, 0, 0
    ref_len = None
    with open(clean_path, 'w') as f_clean, \
         open(log_path, 'w', newline='') as f_log, \
         open(order_path, 'w') as f_order:
        log_writer = csv.writer(f_log)
        log_writer.writerow(['full_header', 'error_reason'])
        for record in SeqIO.parse(input_path, "fasta"):
            processed += 1
            error_reason = None
            seq_id = record.id.split('|')[0]
            if record.id.count('|') < 1:
                error_reason = "Malformed header (expected at least 1 pipe)"
            else:
                if ref_len is None: ref_len = len(record.seq)
                elif len(record.seq) != ref_len:
                    error_reason = f"Incorrect length (is {len(record.seq)}, expected {ref_len})"
            if error_reason:
                errors += 1
                log_writer.writerow([record.id, error_reason])
            else:
                written += 1
                f_clean.write(f">{record.id}\n{record.seq}\n")
                f_order.write(f"{seq_id}\n")
    print("\n--- Verification and Filtering Complete ---")
    print(f"Total sequences processed: {processed}")
    print(f"Sequences written to new file: {written}")
    print(f"Sequences discarded (errors): {errors}")

verify_and_filter_fasta(INPUT_FASTA, CLEAN_FASTA_OUTPUT, LOG_FILE_OUTPUT, ORDER_INDEX_OUTPUT)

# The following would be uncommented in Colab, but in its own Cell
#  of shell code only (not mixed with Python code).
# Commented out IPython magic to ensure Python compatibility.
# %%shell
# RELEASE_DATE="2025/07/05"
# BASE_URL="http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/${RELEASE_DATE}/"
# METADATA_FILE_NAME="public-${RELEASE_DATE//\//-}.metadata.tsv"
# 
# echo "Downloading and decompressing data file from ${RELEASE_DATE}..."
# 
# # Download and rename metadata file
# wget -nc "${BASE_URL}/${METADATA_FILE_NAME}.gz" -O public-latest.metadata.tsv.gz
# gunzip -f public-latest.metadata.tsv.gz

# This cell converts our clean text-based alignment into the high-performance memory-mapped binary format.
# This may also take 10-15 minutes.

# Cell 3: Pre-processing (FASTA to Binary)
import numpy as np
import os
from Bio import SeqIO

print("\n--- Starting Pre-processing: FASTA to Binary ---")

# --- Configuration ---
CLEAN_FASTA_LOCAL = "verified_sequences.fas"
NUMPY_BINARY_FILE = "alignment.mmap"
# --------------------

def fasta_to_numpy_binary(fasta_path, output_path):
    print("First pass to get dimensions...")
    records = list(SeqIO.parse(fasta_path, "fasta"))
    num_sequences = len(records)
    seq_length = len(records[0].seq)
    print(f"Alignment dimensions: {num_sequences} sequences, {seq_length} sites.")

    mmap_array = np.memmap(output_path, dtype=np.uint8, mode='w+', shape=(num_sequences, seq_length))
    nt_map = str.maketrans("ACGTN", "01234")
    print("Second pass to write data...")
    for i, record in enumerate(records):
        processed_sequence = str(record.seq).upper().translate(nt_map)
        mmap_array[i, :] = np.frombuffer(processed_sequence.encode(), dtype=np.uint8) - ord('0')
    mmap_array.flush()
    print(f"Successfully created binary file at: {output_path}")
    return output_path, num_sequences, seq_length

binary_file, n_seq, n_len = fasta_to_numpy_binary(CLEAN_FASTA_LOCAL, NUMPY_BINARY_FILE)

# This portion performs the metadata alignment and then launches the main
# calculation loop. It incorporates both π metrics. 

# Cell 4: Main Analysis - Metadata Alignment and Partitioned Pi Calculation
import numpy as np
import pandas as pd
import gc
import os

print("\n--- Phase 3: Main Analysis ---")

# --- 1. Memory-Efficient Metadata Alignment ---
print("Step 1: Aligning Metadata with Sequence Data...")
ORDER_INDEX_FILE = "verified_sequence_order.txt"
METADATA_TSV = "public-latest.metadata.tsv"
CHUNK_SIZE = 500000

ordered_ids_set = set(pd.read_csv(ORDER_INDEX_FILE, header=None, names=['strain'])['strain'])
cols_to_use = ['strain', 'pangolin_lineage', 'date']
filtered_chunks = []
with pd.read_csv(METADATA_TSV, sep='\t', usecols=cols_to_use, on_bad_lines='skip', chunksize=CHUNK_SIZE) as reader:
    for chunk in reader:
        chunk['strain_norm'] = chunk['strain'].str.split('|').str[0]
        filtered_chunks.append(chunk[chunk['strain_norm'].isin(ordered_ids_set)])
df_meta_filtered = pd.concat(filtered_chunks, ignore_index=True)
del filtered_chunks; gc.collect()
df_meta_filtered.drop_duplicates(subset=['strain_norm'], keep='first', inplace=True)
ordered_ids_list = pd.read_csv(ORDER_INDEX_FILE, header=None, names=['strain_norm'])['strain_norm'].tolist()
df_meta_aligned = df_meta_filtered.set_index('strain_norm').loc[ordered_ids_list].reset_index()
del df_meta_filtered, ordered_ids_list; gc.collect()
print("Metadata alignment complete.")

# --- 2. Define Analytical Partitions ---
print("Step 2: Defining analytical partitions...")
df_meta_aligned['date'] = pd.to_datetime(df_meta_aligned['date'], errors='coerce')
df_meta_aligned.dropna(subset=['date'], inplace=True)
# Define our epochs for the JN.1 sweep analysis
df_meta_aligned['time_bin'] = pd.to_datetime(df_meta_aligned['date']).dt.to_period('W') # Weekly bins for finer resolution
grouped = df_meta_aligned.groupby(['pangolin_lineage', 'time_bin'])
print(f"Created {len(grouped)} groups for analysis.")

# --- 3. Partitioned π Calculation Loop ---
print("Step 3: Starting Partitioned π Calculation Loop (this will take several hours)...")
# This is our enhanced, validated function
def calculate_both_pi(alignment_slice):
    n_seq, n_sites = alignment_slice.shape
    if n_seq < 2: return 0.0, 0.0
    pi_sum_at_valid_sites = 0.0
    sites_compared = 0
    for i in range(n_sites):
        col = alignment_slice[:, i]
        valid_bases = col[col < 4]
        n_k = len(valid_bases)
        if n_k < 2: continue
        sites_compared += 1
        counts = np.bincount(valid_bases, minlength=4)
        sum_freq_sq = np.sum((counts / n_k)**2)
        pi_site = (n_k / (n_k - 1)) * (1 - sum_freq_sq)
        pi_sum_at_valid_sites += pi_site
    pi_valid_sites = pi_sum_at_valid_sites / sites_compared if sites_compared > 0 else 0.0
    pi_total_length = pi_sum_at_valid_sites / n_sites if n_sites > 0 else 0.0
    return pi_valid_sites, pi_total_length

alignment = np.memmap(NUMPY_BINARY_FILE, dtype=np.uint8, mode='r', shape=(n_seq, n_len))
results = []
group_counter = 0
for name, group_df in grouped:
    group_counter += 1
    if len(group_df) < 2: continue
    indices = group_df.index.tolist()
    alignment_slice = alignment[indices, :]
    pi_valid, pi_total = calculate_both_pi(alignment_slice)
    results.append({
        'lineage': name[0],
        'time_bin': str(name[1]),
        'pi_valid_sites': pi_valid,
        'pi_total_length': pi_total,
        'n_samples': len(group_df)
    })
    if group_counter % 500 == 0:
        print(f"  Processed {group_counter}/{len(grouped)} groups...")

print("\n--- Calculation Complete ---")
results_df = pd.DataFrame(results)
print("Sample of results:")
print(results_df.head())

# --- 4. Save Results ---
DRIVE_PATH = "/content/drive/MyDrive/"
RESULTS_FILE = "jn1_nucleotide_diversity_results.csv"
results_df.to_csv(RESULTS_FILE, index=False)
!cp "{RESULTS_FILE}" "{DRIVE_PATH}"
print(f"\nResults saved to {RESULTS_FILE} and copied to Google Drive.")

"""END OF PART 2 AND 3"""

# Cell for Phase 4 - Analysis and Plotting Setup

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from google.colab import drive
import os

# Mount Google Drive
drive.mount('/content/drive')

print("--- Phase 4: Analysis and Visualization ---")

# --- 1. Load and Filter Data ---
print("Step 1: Loading and filtering results...")
DRIVE_PATH = "/content/drive/MyDrive/"
RESULTS_FILE = "jn1_nucleotide_diversity_results.csv"

try:
    df = pd.read_csv(f"{DRIVE_PATH}{RESULTS_FILE}")
    print(f"Successfully loaded {len(df)} raw data points.")
except FileNotFoundError:
    print(f"FATAL: Results file not found at {DRIVE_PATH}{RESULTS_FILE}")
    # Stop execution if file not found

# Convert time_bin string to a proper datetime object for plotting
df['time_bin_start'] = pd.to_datetime(df['time_bin'].str.split('/').str[0])

# Apply a filter for statistical robustness
MIN_SAMPLES = 10
df_filtered = df[df['n_samples'] >= MIN_SAMPLES].copy()
print(f"Filtered data to {len(df_filtered)} points (n_samples >= {MIN_SAMPLES}).")

# --- 2. Calculate Global Nucleotide Diversity ---
print("Step 2: Calculating global nucleotide diversity per time bin...")

def weighted_average(group, avg_name, weight_name):
    """Helper function to calculate weighted average."""
    d = group[avg_name]
    w = group[weight_name]
    return (d * w).sum() / w.sum()

# Group by time bin and calculate the weighted average pi for each week
global_pi = df_filtered.groupby('time_bin_start').apply(
    weighted_average, 'pi_valid_sites', 'n_samples'
).to_frame(name='global_pi_valid')

global_pi['global_pi_total'] = df_filtered.groupby('time_bin_start').apply(
    weighted_average, 'pi_total_length', 'n_samples'
)

print("Global diversity calculation complete. Ready for plotting.")
print("\nSample of global diversity data:")
print(global_pi.head())

# --- Ready for Plotting ---
# The 'global_pi' DataFrame is now ready to be used for our primary visualizations.

# Cell for Phase 4 - Primary Visualization of the JN.1 Sweep

# --- 1. Plotting Setup ---
sns.set_style("whitegrid")
plt.figure(figsize=(16, 8))

# --- 2. Focus on the Relevant Time Period ---
# We slice the data to start from mid-2023 to focus on the JN.1 event.
start_date = '2023-07-01'
global_pi_slice = global_pi[global_pi.index >= start_date]

# --- 3. Create the Time-Series Plot ---
plt.plot(global_pi_slice.index, global_pi_slice['global_pi_valid'],
         marker='o', linestyle='-', markersize=4, label='Global Nucleotide Diversity (π)')

# --- 4. Annotate the Plot for Clarity ---
# Add vertical lines to demarcate our defined epochs
plt.axvline(pd.to_datetime('2023-11-01'), color='r', linestyle='--', alpha=0.7, label='Approx. Start of JN.1 Sweep')
plt.axvline(pd.to_datetime('2024-02-01'), color='g', linestyle='--', alpha=0.7, label='Approx. End of Sweep (JN.1 Dominant)')

# Add labels and title
plt.title('Global SARS-CoV-2 Nucleotide Diversity During the JN.1 Sweep', fontsize=18, pad=20)
plt.xlabel('Date', fontsize=12)
plt.ylabel('Global Nucleotide Diversity (π)', fontsize=12)
plt.legend()
plt.ylim(bottom=0) # Ensure the y-axis starts at 0

# Improve date formatting on the x-axis
plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%b %Y'))
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()

# Cell for Phase 4 - Comparative Plot of Both Pi Metrics

# --- 1. Plotting Setup ---
sns.set_style("whitegrid")
plt.figure(figsize=(16, 8))

# --- 2. Focus on the Relevant Time Period ---
start_date = '2023-07-01'
global_pi_slice = global_pi[global_pi.index >= start_date]

# --- 3. Create the Comparative Time-Series Plot ---
# Plot our primary, validated metric
plt.plot(global_pi_slice.index, global_pi_slice['global_pi_valid'],
         marker='o', linestyle='-', markersize=4, label='π (normalized by valid sites)')

# Plot the alternate metric for comparison
plt.plot(global_pi_slice.index, global_pi_slice['global_pi_total'],
         marker='x', linestyle='--', markersize=4, alpha=0.8, label='π (normalized by total length)')

# --- 4. Annotate the Plot for Clarity ---
plt.axvline(pd.to_datetime('2023-11-01'), color='r', linestyle='--', alpha=0.7, label='Approx. Start of JN.1 Sweep')
plt.axvline(pd.to_datetime('2024-02-01'), color='g', linestyle='--', alpha=0.7, label='Approx. End of Sweep')

# Add labels and title
plt.title('Comparison of Global SARS-CoV-2 Nucleotide Diversity Metrics During the JN.1 Sweep', fontsize=18, pad=20)
plt.xlabel('Date', fontsize=12)
plt.ylabel('Global Nucleotide Diversity (π)', fontsize=12)
plt.legend()
plt.ylim(bottom=0)

plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%b %Y'))
plt.xticks(rotation=45, ha='right')

plt.tight_layout()
plt.show()

# Cell for Phase 4 - Intra-Lineage Diversity Plot for JN.1 (Definitive and Final)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# --- 1. Plotting Setup ---
sns.set_style("whitegrid")
plt.figure(figsize=(16, 8))

# --- 2. Filter for JN.1 and its sublineages using the correct notation ---
# The lineage column in the metadata uses '.', not '_'.
df_jn1_only = df_filtered[df_filtered['lineage'].str.startswith('JN.')].copy()
# ---------------------------------------------

# Check if the filter produced any data
if df_jn1_only.empty:
    print("Warning: No data found for JN.1 lineages after filtering. Cannot generate plot.")
else:
    # --- 3. Calculate Weighted Average ---
    # To calculate the weighted average robustly, we first calculate the numerator part for each row.
    df_jn1_only['pi_x_samples'] = df_jn1_only['pi_valid_sites'] * df_jn1_only['n_samples']

    # Now, we group by the time bin and sum the necessary components.
    grouped_sums = df_jn1_only.groupby('time_bin_start').agg({
        'pi_x_samples': 'sum',
        'n_samples': 'sum'
    })

    # The weighted average is the sum of (pi * n) divided by the sum of n.
    jn1_intra_pi = (grouped_sums['pi_x_samples'] / grouped_sums['n_samples']).to_frame(name='intra_jn1_pi')

    # --- 4. Create the Time-Series Plot ---
    plt.plot(jn1_intra_pi.index, jn1_intra_pi['intra_jn1_pi'],
             marker='o', linestyle='-', markersize=4, label='Intra-Lineage Diversity (π) within JN.1 Family')

    # --- 5. Annotate the Plot for Clarity ---
    plt.title('Post-Sweep Diversification within the JN.1 Lineage', fontsize=18, pad=20)
    plt.xlabel('Date', fontsize=12)
    plt.ylabel('Nucleotide Diversity (π)', fontsize=12)
    plt.legend()
    plt.ylim(bottom=0)

    # Focus on the relevant time period
    if not jn1_intra_pi.empty:
        plt.xlim(pd.to_datetime('2023-10-01'), jn1_intra_pi.index.max() + pd.Timedelta(days=30))

    plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%b %Y'))
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    plt.show()
