# analysis_pipeline.py
#
# Description:
# This script performs a comprehensive genomic analysis of SARS-CoV-2 lineages.
# It is designed to reproduce the findings presented in two reports:
#   1. A comparative nucleotide diversity (Pi) analysis between JN.1 and BA.2.86.
#   2. A temporal analysis of key selective sweeps (L455S and F456L) within the JN.1 lineage.
#
# The pipeline proceeds in the following stages:
#   1. Initial Setup: Verifies a master alignment and converts it to a binary format.
#   2. Data Loading: Loads metadata, GFF annotations, and the binary alignment.
#   3. Analysis for Report 1 (Comparative Diversity):
#      - Calculates sliding window nucleotide diversity for JN.1 and BA.2.86.
#      - Generates a comparative plot of the diversity profiles.
#      - Performs a codon-level analysis to identify key differentiating mutations.
#   4. Analysis for Report 2 (Temporal Sweeps):
#      - Calculates and plots the temporal frequency of the L455S mutation.
#      - Calculates and plots the temporal frequency of the F456L mutation.
#
# Prerequisites:
# - A `public-latest.metadata.tsv` file from UShER.
# - An `all_clades_aligned.fas` file (or an archive containing it).
# - A `NC_045512_2.gff3` annotation file.
# - A `NC_045512v2.fa` reference genome file.

import os
import re
import csv
import gc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

# --- 0. Master Configuration ---
print("--- Initializing Master Configuration ---")
# File Paths
METADATA_TSV = "public-latest.metadata.tsv"
RAW_FASTA = "all_clades_aligned.fas"
CLEAN_FASTA = "verified_sequences.fas"
ORDER_INDEX_FILE = "verified_sequence_order.txt"
NUMPY_BINARY_FILE = "alignment.mmap"
GFF_FILE = "NC_045512_2.gff3"
REF_FASTA_FILE = "NC_045512v2.fa"

# Analysis Parameters
CHUNK_SIZE = 500000  # For reading large metadata files
WINDOW_SIZE = 500
STEP_SIZE = 50


# --- 1. Utility Functions ---
def verify_and_filter_fasta(input_path, clean_path, order_path):
    """Verifies sequence length and header format, writing valid sequences to a new file."""
    print(f"\n--- Verifying FASTA file: {input_path} ---")
    if not os.path.exists(input_path):
        print(f"FATAL: Input FASTA not found at {input_path}. Skipping verification.")
        return False
    processed, written, errors = 0, 0, 0
    ref_len = None
    with open(clean_path, 'w') as f_clean, open(order_path, 'w') as f_order:
        for record in SeqIO.parse(input_path, "fasta"):
            processed += 1
            error_reason = None
            seq_id = record.id.split('|')[0]
            if record.id.count('|') < 1: error_reason = "Malformed header"
            else:
                if ref_len is None: ref_len = len(record.seq)
                elif len(record.seq) != ref_len: error_reason = "Incorrect length"
            if error_reason:
                errors += 1
            else:
                written += 1
                f_clean.write(f">{record.id}\n{record.seq}\n")
                f_order.write(f"{seq_id}\n")
    print(f"Verification Complete: {written} sequences written, {errors} discarded.")
    return True

def fasta_to_numpy_binary(fasta_path, output_path):
    """Converts a FASTA file to a memory-mapped NumPy binary array."""
    print(f"\n--- Converting FASTA to Binary: {fasta_path} ---")
    if not os.path.exists(fasta_path):
        print(f"FATAL: Clean FASTA not found at {fasta_path}. Skipping binary conversion.")
        return None, 0, 0
    records = list(SeqIO.parse(fasta_path, "fasta"))
    num_sequences, seq_length = len(records), len(records[0].seq)
    print(f"Alignment dimensions: {num_sequences} sequences, {seq_length} sites.")
    mmap_array = np.memmap(output_path, dtype=np.uint8, mode='w+', shape=(num_sequences, seq_length))
    nt_map = str.maketrans("ACGTN", "01234")
    for i, record in enumerate(records):
        processed_seq = str(record.seq).upper().translate(nt_map)
        mmap_array[i, :] = np.frombuffer(processed_seq.encode(), dtype=np.uint8) - ord('0')
    mmap_array.flush()
    print(f"Successfully created binary file at: {output_path}")
    return output_path, num_sequences, seq_length

def parse_gff3_to_df(gff_file_path):
    """Parses a GFF3 file and returns a DataFrame with gene annotations."""
    gff_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    try:
        df = pd.read_csv(gff_file_path, sep='\t', comment='#', header=None, names=gff_columns)
    except FileNotFoundError: return None
    df_genes = df[df['type'] == 'gene'].copy()
    df_genes['gene_name'] = df_genes['attributes'].apply(lambda x: re.search(r'Name=([^;]+)', x).group(1) if re.search(r'Name=([^;]+)', x) else None)
    return df_genes[['gene_name', 'start', 'end']].dropna().reset_index(drop=True)

def load_and_align_metadata(metadata_path, order_path):
    """Loads and filters metadata to match the order of sequences in the alignment."""
    print("\n--- Loading and Aligning Metadata ---")
    ordered_ids_set = set(pd.read_csv(order_path, header=None, names=['strain'])['strain'])
    cols_to_use = ['strain', 'pangolin_lineage', 'date', 'country']
    filtered_chunks = []
    with pd.read_csv(metadata_path, sep='\t', usecols=cols_to_use, on_bad_lines='skip', chunksize=CHUNK_SIZE) as reader:
        for chunk in reader:
            chunk['strain_norm'] = chunk['strain'].str.split('|').str[0]
            filtered_chunks.append(chunk[chunk['strain_norm'].isin(ordered_ids_set)])
    df_meta_filtered = pd.concat(filtered_chunks, ignore_index=True)
    df_meta_filtered.drop_duplicates(subset=['strain_norm'], keep='first', inplace=True)
    ordered_ids_list = pd.read_csv(order_path, header=None, names=['strain_norm'])['strain_norm'].tolist()
    df_meta_aligned = df_meta_filtered.set_index('strain_norm').loc[ordered_ids_list].reset_index()
    print("Metadata alignment complete.")
    return df_meta_aligned


# --- 2. Main Execution Block ---
if __name__ == "__main__":
    # --- STAGE 1: SETUP AND PRE-PROCESSING ---
    # These steps only need to be run once to prepare the data.
    if not os.path.exists(CLEAN_FASTA):
        verify_and_filter_fasta(RAW_FASTA, CLEAN_FASTA, ORDER_INDEX_FILE)
    if not os.path.exists(NUMPY_BINARY_FILE):
        fasta_to_numpy_binary(CLEAN_FASTA, NUMPY_BINARY_FILE)

    # --- STAGE 2: LOAD ALL NECESSARY DATA FOR ANALYSIS ---
    print("\n--- Loading All Analytical Data ---")
    df_meta_aligned = load_and_align_metadata(METADATA_TSV, ORDER_INDEX_FILE)
    gene_annotations_df = parse_gff3_to_df(GFF_FILE)
    ref_record = SeqIO.read(REF_FASTA_FILE, "fasta")
    ref_seq = ref_record.seq

    # Load memory-mapped alignment
    n_seq = len(df_meta_aligned)
    file_size = os.path.getsize(NUMPY_BINARY_FILE)
    n_len = file_size // n_seq
    alignment = np.memmap(NUMPY_BINARY_FILE, dtype=np.uint8, mode='r', shape=(n_seq, n_len))
    nt_map_rev = {0: 'A', 1: 'C', 2: 'G', 3: 'T', 4: 'N'}
    print(f"Successfully loaded memory-mapped alignment with shape: {alignment.shape}")

    # --- STAGE 3: ANALYSIS FOR REPORT 1 (COMPARATIVE DIVERSITY) ---
    print("\n\n--- STARTING ANALYSIS FOR REPORT 1: COMPARATIVE DIVERSITY ---")

    # Define lineage groups
    is_jn1 = df_meta_aligned['pangolin_lineage'].str.startswith('JN.1', na=False)
    is_ba286_family = df_meta_aligned['pangolin_lineage'].str.startswith('BA.2.86', na=False)
    jn1_indices = df_meta_aligned[is_jn1].index.tolist()
    ba286_parent_indices = df_meta_aligned[is_ba286_family & ~is_jn1].index.tolist()
    print(f"Found {len(jn1_indices)} JN.1 sequences and {len(ba286_parent_indices)} BA.2.86 (Parent) sequences.")

    # A. Sliding Window Nucleotide Diversity (Pi)
    def calculate_pi_for_window(alignment_slice):
        n_s, n_sites = alignment_slice.shape
        if n_s < 2: return 0.0
        pi_sum = 0.0
        sites_compared = 0
        for i in range(n_sites):
            col = alignment_slice[:, i]
            valid_bases = col[col < 4]
            n_k = len(valid_bases)
            if n_k < 2: continue
            sites_compared += 1
            counts = np.bincount(valid_bases, minlength=4)
            pi_sum += (n_k / (n_k - 1)) * (1 - np.sum((counts / n_k)**2))
        return pi_sum / sites_compared if sites_compared > 0 else 0.0

    def run_sliding_window_analysis(indices, lineage_name):
        print(f"\nCalculating Pi for {lineage_name}...")
        results = []
        lineage_alignment = alignment[indices, :]
        for start in range(0, n_len - WINDOW_SIZE + 1, STEP_SIZE):
            window = lineage_alignment[:, start:start+WINDOW_SIZE]
            pi_value = calculate_pi_for_window(window)
            results.append({'window_midpoint': start + (WINDOW_SIZE // 2), 'pi': pi_value})
        df = pd.DataFrame(results)
        df['lineage'] = lineage_name
        return df

    df_pi_jn1 = run_sliding_window_analysis(jn1_indices, "JN.1")
    df_pi_ba286 = run_sliding_window_analysis(ba286_parent_indices, "BA.2.86 (Parent)")
    combined_pi_df = pd.concat([df_pi_jn1, df_pi_ba286], ignore_index=True)

    # B. Plot Comparative Diversity
    print("\nGenerating comparative diversity plot...")
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.figure(figsize=(20, 8))
    ax = sns.lineplot(data=combined_pi_df, x='window_midpoint', y='pi', hue='lineage', linewidth=1.5, zorder=2)
    y_max = combined_pi_df['pi'].max()
    for _, row in gene_annotations_df.iterrows():
        plt.axvspan(row['start'], row['end'], color='grey', alpha=0.08, zorder=1)
        plt.text(x=(row['start'] + row['end']) / 2, y=y_max * 1.05, s=row['gene_name'],
                 ha='center', va='bottom', fontsize=9, rotation=90)
    plt.title('Comparative Analysis of Nucleotide Diversity: JN.1 vs. BA.2.86', fontsize=16, pad=20)
    plt.xlabel('Genomic Position (Window Midpoint)', fontsize=12)
    plt.ylabel('Nucleotide Diversity (π)', fontsize=12)
    plt.xlim(0, 30000); plt.ylim(bottom=0, top=y_max * 1.2)
    plt.legend(title='Lineage', loc='upper left')
    plt.tight_layout()
    plt.savefig("jn1_ba286_sliding_window.png")
    plt.show()

    # C. Codon-level analysis for key mutations
    print("\nPerforming codon-level analysis for key genes...")
    mutation_results = []
    genes_of_interest = ['S', 'ORF3a', 'N']
    for _, gene_row in gene_annotations_df[gene_annotations_df['gene_name'].isin(genes_of_interest)].iterrows():
        gene_name, start, end = gene_row['gene_name'], gene_row['start'] - 1, gene_row['end']
        for i in range(0, (end - start) - 2, 3):
            codon_num = (i // 3) + 1
            ref_aa = str(ref_seq[start+i:start+i+3].translate())
            
            def get_consensus_aa(indices):
                counts = Counter()
                for sample_idx in indices:
                    codon_n = alignment[sample_idx, start+i:start+i+3]
                    if 4 in codon_n: continue
                    aa = str(Seq("".join([nt_map_rev[n] for n in codon_n])).translate())
                    counts[aa] += 1
                if not counts: return None, 0.0
                consensus = counts.most_common(1)[0]
                return consensus[0], consensus[1] / sum(counts.values())

            jn1_aa, jn1_freq = get_consensus_aa(jn1_indices)
            ba286_aa, ba286_freq = get_consensus_aa(ba286_parent_indices)

            if jn1_aa and jn1_aa != ref_aa and (jn1_aa != ba286_aa or jn1_freq > ba286_freq + 0.1):
                mutation_results.append({
                    'Gene': gene_name, 'Mutation': f"{ref_aa}{codon_num}{jn1_aa}",
                    'JN.1 Freq': round(jn1_freq, 3), 'BA.2.86 Freq': round(ba286_freq, 3)
                })
    print("Identified non-synonymous mutations enriched in JN.1:")
    print(pd.DataFrame(mutation_results).to_string())


    # --- STAGE 4: ANALYSIS FOR REPORT 2 (TEMPORAL SWEEPS) ---
    print("\n\n--- STARTING ANALYSIS FOR REPORT 2: TEMPORAL SWEEPS ---")
    
    # Prepare temporal data for JN.1
    temporal_jn1_meta = df_meta_aligned.loc[jn1_indices].copy()
    temporal_jn1_meta['date'] = pd.to_datetime(temporal_jn1_meta['date'], errors='coerce')
    temporal_jn1_meta.dropna(subset=['date'], inplace=True)
    temporal_jn1_meta['week'] = temporal_jn1_meta['date'].dt.to_period('W').apply(lambda p: p.start_time)

    def run_temporal_analysis(gene_name, codon_number, target_aa, title, filename):
        print(f"\nCalculating temporal frequency for {gene_name}:{target_aa}{codon_number}...")
        time_series_data = []
        try:
            gene_info = gene_annotations_df[gene_annotations_df['gene_name'] == gene_name].iloc[0]
            codon_start_abs = (gene_info['start'] - 1) + ((codon_number - 1) * 3)
            
            for week, group in temporal_jn1_meta.groupby('week'):
                indices = group.index.tolist()
                codons_numeric = alignment[indices, codon_start_abs:codon_start_abs+3]
                
                target_count, total_valid = 0, 0
                for codon_num in codons_numeric:
                    if 4 in codon_num: continue
                    total_valid += 1
                    aa = str(Seq("".join([nt_map_rev[n] for n in codon_num])).translate())
                    if aa == target_aa: target_count += 1
                
                if total_valid > 0:
                    time_series_data.append({'week': week, 'frequency': target_count / total_valid})
            
            freq_df = pd.DataFrame(time_series_data)
            
            # Plotting
            plt.figure(figsize=(14, 7))
            sns.lineplot(data=freq_df, x='week', y='frequency', marker='o', linestyle='-')
            plt.title(title, fontsize=16, pad=20)
            plt.xlabel('Collection Date (Week)', fontsize=12)
            plt.ylabel(f'Frequency of {target_aa} at Spike Position {codon_number}', fontsize=12)
            plt.ylim(0, 1.05)
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(filename)
            plt.show()
        except (IndexError, KeyError) as e:
            print(f"FATAL: Could not perform temporal analysis. Reason: {e}")

    # A. Plot for lineage-defining L455S mutation
    run_temporal_analysis(
        gene_name='S', codon_number=455, target_aa='S',
        title='Frequency of Lineage-Defining S455 Allele in JN.1',
        filename='jn1_s455_frequency.png'
    )

    # B. Plot for emerging F456L mutation
    run_temporal_analysis(
        gene_name='S', codon_number=456, target_aa='L',
        title='Frequency of Emerging F456L Allele in the JN.1 Background',
        filename='f456l_frequency.png'
    )

    print("\n\n--- All Analyses Complete ---")