import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import Counter
from itertools import product
from scipy.stats import pearsonr
from scipy.signal import find_peaks
import shutil

# Step 1: Read the genome from a FASTA file
def read_genome(filename):
    """Reads a single genome sequence from a FASTA file."""
    record = next(SeqIO.parse(filename, "fasta"))
    return str(record.seq)

# Step 2: Calculate TNF (Tetranucleotide Frequency) for a given sequence
def calculate_tnf(sequence):
    """Calculates the frequency of all 4-mers in a DNA sequence."""
    tnf = Counter()
    for i in range(len(sequence) - 3):
        kmer = sequence[i:i + 4].upper()
        if all(base in "ATGC" for base in kmer):
            tnf[kmer] += 1
    return tnf

# Step 3: Normalize TNF to create consistent k-mer vectors
def normalize_tnf(tnf, all_kmers):
    """Normalizes TNF counts to create a frequency vector."""
    normalized = {kmer: tnf.get(kmer, 0) for kmer in all_kmers}
    total = sum(normalized.values())
    return {k: (v / total) if total > 0 else 0 for k, v in normalized.items()}

# Step 4: Calculate Pearson correlations for each chunk size
def calculate_chunk_correlations(genome_seq, chunk_size):
    """Divides the genome into chunks and calculates the Pearson correlation of each chunk's TNF against the whole genome's TNF."""
    all_kmers = [''.join(p) for p in product('ATGC', repeat=4)]
    genome_tnf = normalize_tnf(calculate_tnf(genome_seq), all_kmers)
    genome_tnf_vector = np.array(list(genome_tnf.values()), dtype=np.float64)

    correlations = []
    for i in range(0, len(genome_seq), chunk_size):
        chunk_seq = genome_seq[i:i + chunk_size]
        chunk_tnf = normalize_tnf(calculate_tnf(chunk_seq), all_kmers)
        chunk_tnf_vector = np.array(list(chunk_tnf.values()), dtype=np.float64)

        # Add a small epsilon to avoid zero variance issues
        genome_tnf_vector += 1e-9
        chunk_tnf_vector += 1e-9

        correlation, _ = pearsonr(genome_tnf_vector, chunk_tnf_vector)
        correlations.append(correlation)

    return correlations

# Step 5: Calculate cumulative sum (cumsum)
def calculate_cumsum(correlations):
    """Calculates the cumulative sum of the correlation values."""
    return np.cumsum(correlations)

# Step 6: Calculate the first derivative
def calculate_first_derivative(cumsum_values):
    """Calculates the first derivative of the cumsum plot to identify rate of change."""
    return np.gradient(cumsum_values) if len(cumsum_values) >= 2 else np.array([])

# Step 7: Calculate sliding window variance
def sliding_window_variance(data, window_size=5):
    """Calculates the variance in a sliding window across the data to find unstable regions."""
    return [np.var(data[i:i + window_size]) for i in range(len(data) - window_size + 1)]

# Step 8: Detect genomic islands with adaptive thresholds
def detect_genomic_islands(variance_data, chunk_size):
    """Detects peaks in variance data using adaptive thresholds to identify potential GIs."""
    if len(variance_data) == 0:
        print("Warning: Variance data is empty; skipping genomic island detection.")
        return [], []

    # Adaptive thresholds based on chunk size
    if chunk_size == 5000:
        primary_threshold = np.percentile(variance_data, 90)
        secondary_threshold = np.mean(variance_data) + np.std(variance_data)
    else:  # chunk_size == 2000
        primary_threshold = np.percentile(variance_data, 80)
        secondary_threshold = np.mean(variance_data) + (0.5 * np.std(variance_data))

    primary_peaks, _ = find_peaks(variance_data, height=primary_threshold)
    secondary_peaks, _ = find_peaks(variance_data, height=secondary_threshold)
    
    all_peaks = np.unique(np.concatenate((primary_peaks, secondary_peaks)))
    peak_types = ['Primary' if peak in primary_peaks else 'Secondary' for peak in all_peaks]

    return all_peaks, peak_types

# Step 9: Extend genomic island boundaries
def extend_island_boundaries(variance_data, peaks):
    """Extends the boundaries of a GI from its peak outwards until the variance returns to baseline."""
    baseline = np.mean(variance_data)
    extended_islands = []
    for peak in peaks:
        start = peak
        while start > 0 and variance_data[start] > baseline:
            start -= 1
        end = peak
        while end < len(variance_data) - 1 and variance_data[end] > baseline:
            end += 1
        extended_islands.append((start, end))
    return extended_islands

# Step 10: Merge overlapping or adjacent genomic islands
def merge_overlapping_islands(islands):
    """Merges GI regions that overlap or are adjacent to each other."""
    if not islands:
        return []
    islands = sorted(islands, key=lambda x: x[0])
    merged_islands = [islands[0]]
    for current_start, current_end in islands[1:]:
        last_start, last_end = merged_islands[-1]
        if current_start <= last_end + 1:
            merged_islands[-1] = (last_start, max(last_end, current_end))
        else:
            merged_islands.append((current_start, current_end))
    return merged_islands

# NEW Step 11: Find Direct Repeats flanking a Genomic Island
def find_direct_repeats_for_island(genome_seq, island_start_bp, island_end_bp,
                                     search_window=500, min_repeat_len=10, max_repeat_len=50):
    """
    Searches for direct repeats in the flanking regions of a genomic island.

    Returns:
        A tuple (bool, dict): 
        - True if a repeat is found, False otherwise.
        - A dictionary with details of the found repeat, or None.
    """
    # Define upstream and downstream search regions, handling genome boundaries
    upstream_start = max(0, island_start_bp - search_window)
    upstream_region = genome_seq[upstream_start:island_start_bp]
    
    downstream_end = min(len(genome_seq), island_end_bp + search_window)
    downstream_region = genome_seq[island_end_bp:downstream_end]

    if not upstream_region or not downstream_region:
        return False, None

    # Search for repeats, starting with longer ones for higher significance
    for length in range(max_repeat_len, min_repeat_len - 1, -1):
        for i in range(len(upstream_region) - length + 1):
            repeat_candidate = upstream_region[i:i + length]
            downstream_pos = downstream_region.find(repeat_candidate)
            
            if downstream_pos != -1:
                # Repeat found!
                repeat_details = {
                    "seq": repeat_candidate,
                    "length": length,
                    "upstream_pos": upstream_start + i,
                    "downstream_pos": island_end_bp + downstream_pos
                }
                return True, repeat_details
                
    return False, None

# Step 12: Save genomic islands to FASTA files
def save_islands_to_fasta_from_csv(genome_seq, csv_file, sequence_id, output_dir, chunk_size_kb):
    """Saves the final identified island sequences to a FASTA file."""
    if not os.path.exists(csv_file):
        print(f"Warning: CSV file not found, skipping FASTA generation: {csv_file}")
        return
        
    df = pd.read_csv(csv_file)
    fasta_records = []
    for index, row in df.iterrows():
        nucleotide_start = int(row['Start'])
        nucleotide_end = int(row['End'])
        island_type = row.get('Peak Type', 'Island')
        island_seq = genome_seq[nucleotide_start:nucleotide_end]
        header = f"{sequence_id}_start{nucleotide_start}_end{nucleotide_end}_{island_type}"
        fasta_records.append(SeqRecord(Seq(island_seq), id=header, description=""))

    fasta_filename = os.path.join(output_dir, f"{sequence_id}_identified_islands_{chunk_size_kb}kb.fasta")
    SeqIO.write(fasta_records, fasta_filename, "fasta")
    print(f"FASTA file for {chunk_size_kb}kb islands saved at {fasta_filename}")

# MODIFIED Step 13: Plot results and mark islands with confidence levels
def plot_results_with_extended_islands(correlations, cumsum_values, first_derivative, variances, extended_islands, peak_types, output_dir, chunk_size):
    """Plots all analysis steps and highlights the final GIs with different colors based on confidence."""
    plt.figure(figsize=(16, 12))
    
    # Plot 1: Pearson Correlations
    plt.subplot(4, 1, 1)
    plt.plot(correlations, marker='.', linestyle='-', markersize=4, label=f'Pearson Correlations ({chunk_size}bp Chunks)')
    plt.title(f'Analysis for {os.path.basename(output_dir)} ({chunk_size}bp Chunks)', fontsize=16)
    plt.ylabel('Correlation')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()

    # Plot 2: Cumulative Sum
    plt.subplot(4, 1, 2)
    plt.plot(cumsum_values, marker='.', linestyle='-', markersize=4, color='g', label='Cumulative Sum (Cumsum)')
    plt.ylabel('Cumsum')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()

    # Plot 3: First Derivative
    plt.subplot(4, 1, 3)
    if len(first_derivative) > 0:
        plt.plot(first_derivative, marker='.', linestyle='-', markersize=4, color='purple', label='First Derivative')
        plt.ylabel('Rate of Change')
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.legend()

    # Plot 4: Variance with Highlighted GIs
    plt.subplot(4, 1, 4)
    plt.plot(variances, marker='.', linestyle='-', markersize=4, color='black', label='Sliding Window Variance')
    plt.xlabel('Chunk Index', fontsize=12)
    plt.ylabel('Variance')
    plt.grid(True, linestyle='--', alpha=0.6)

    # Use a set to add legend entries only once
    added_labels = set()
    for (start, end), peak_type in zip(extended_islands, peak_types):
        if peak_type == 'High-Confidence':
            color = 'blue'
            label = 'High-Confidence Island (Repeats)'
        elif peak_type == 'Primary':
            color = 'red'
            label = 'Primary Island'
        else:  # Secondary
            color = 'orange'
            label = 'Secondary Island'
        
        if label not in added_labels:
            plt.axvspan(start, end, color=color, alpha=0.4, label=label)
            added_labels.add(label)
        else:
            plt.axvspan(start, end, color=color, alpha=0.4)

    plt.legend()
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    output_path = os.path.join(output_dir, f"genome_analysis_{chunk_size}_islands.png")
    plt.savefig(output_path)
    plt.close()

# MODIFIED Step 14: Output coordinates and repeat info to a CSV
def output_genomic_island_coordinates(island_info_list, chunk_size, output_dir):
    """Outputs the final island coordinates and confidence levels to a CSV file."""
    if not island_info_list:
        print(f"No islands detected for chunk size {chunk_size}, skipping CSV output.")
        return

    df_data = []
    for info in island_info_list:
        details = info['repeat_details']
        df_data.append({
            'Start': info['start_bp'],
            'End': info['end_bp'],
            'Peak Type': info['type'],
            'Direct Repeat Found': 'Yes' if details else 'No',
            'Repeat Sequence': details['seq'] if details else 'N/A',
            'Repeat Length': details['length'] if details else 'N/A',
            'Upstream Position': details['upstream_pos'] if details else 'N/A',
            'Downstream Position': details['downstream_pos'] if details else 'N/A'
        })
    
    df = pd.DataFrame(df_data)
    output_path = os.path.join(output_dir, f"genomic_islands_coordinates_{chunk_size}.csv")
    df.to_csv(output_path, index=False)
    print(f"Genomic islands coordinates for {chunk_size}bp chunks saved to {output_path}")


# MODIFIED Step 15: Process all FASTA records from FASTA file(s)
def process_fasta_files(fasta_paths, base_output_dir, min_length, summary_file="genome_summary.csv"):
    """Main processing loop for FASTA files, including multi-FASTA files."""
    summary_file = os.path.join(base_output_dir, summary_file)
    summary_columns = [
        "Genome ID", "FASTA Header", "Total Genome Size",
        "Number of High-Confidence Islands", "Total High-Confidence Island Size",
        "Number of Primary Islands", "Total Primary Island Size",
        "Number of Secondary Islands", "Total Secondary Island Size"
    ]
    if not os.path.exists(summary_file):
        pd.DataFrame(columns=summary_columns).to_csv(summary_file, index=False)

    chunk_sizes = [5000, 2000]
    used_sequence_ids = set()
    for file_path in fasta_paths:
        fasta_file = os.path.basename(file_path)
        records_found = False
        for record_number, record in enumerate(SeqIO.parse(file_path, "fasta"), start=1):
            records_found = True
            sequence_id = sanitize_sequence_id(record.id)
            if not sequence_id:
                sequence_id = f"{os.path.splitext(fasta_file)[0]}_record_{record_number}"
            original_sequence_id = sequence_id
            duplicate_count = 2
            while sequence_id in used_sequence_ids or os.path.exists(os.path.join(base_output_dir, sequence_id)):
                sequence_id = f"{original_sequence_id}_{duplicate_count}"
                duplicate_count += 1
            used_sequence_ids.add(sequence_id)

            print(f"\n--- Processing Genome: {sequence_id} ---")

            try:
                output_dir = create_output_dir(base_output_dir, sequence_id)
                genome_seq = str(record.seq)
                fasta_header = record.description
                total_genome_size = len(genome_seq)

                if total_genome_size < min_length:
                    print(f"Skipping {sequence_id} (Length: {total_genome_size}bp) - below minimum length {min_length}bp.")
                    continue

                summary_data = {
                    "Genome ID": sequence_id, "FASTA Header": fasta_header, "Total Genome Size": total_genome_size,
                    "Number of High-Confidence Islands": 0, "Total High-Confidence Island Size": 0,
                    "Number of Primary Islands": 0, "Total Primary Island Size": 0,
                    "Number of Secondary Islands": 0, "Total Secondary Island Size": 0
                }

                for chunk_size in chunk_sizes:
                    print(f"Analyzing with {chunk_size}bp chunks...")
                    correlations = calculate_chunk_correlations(genome_seq, chunk_size)
                    cumsum_values = calculate_cumsum(correlations)
                    first_derivative = calculate_first_derivative(cumsum_values)
                    variances = sliding_window_variance(first_derivative, window_size=5)

                    if not variances:
                        print(f"Could not generate variance data for {chunk_size}bp chunks. Skipping.")
                        continue

                    genomic_islands, peak_types = detect_genomic_islands(variances, chunk_size)
                    extended_islands = extend_island_boundaries(variances, genomic_islands)
                    merged_islands = merge_overlapping_islands(extended_islands)

                    final_peak_types = []
                    for start_chunk, end_chunk in merged_islands:
                        original_peak_indices = [p for p in genomic_islands if start_chunk <= p <= end_chunk]
                        original_types = ['Primary' if p in find_peaks(variances, height=np.percentile(variances, 90 if chunk_size==5000 else 80))[0] else 'Secondary' for p in original_peak_indices]
                        final_peak_types.append('Primary' if 'Primary' in original_types else 'Secondary')

                    final_island_info = []
                    for (start_chunk, end_chunk), peak_type in zip(merged_islands, final_peak_types):
                        island_start_bp = start_chunk * chunk_size
                        island_end_bp = min((end_chunk + 1) * chunk_size, total_genome_size)

                        has_repeats, repeat_details = find_direct_repeats_for_island(genome_seq, island_start_bp, island_end_bp)
                        current_type = 'High-Confidence' if has_repeats else peak_type

                        final_island_info.append({
                            "start_chunk": start_chunk, "end_chunk": end_chunk,
                            "start_bp": island_start_bp, "end_bp": island_end_bp,
                            "type": current_type, "repeat_details": repeat_details
                        })

                    plot_islands = [(info['start_chunk'], info['end_chunk']) for info in final_island_info]
                    plot_types = [info['type'] for info in final_island_info]

                    plot_results_with_extended_islands(correlations, cumsum_values, first_derivative, variances, plot_islands, plot_types, output_dir, chunk_size)
                    output_genomic_island_coordinates(final_island_info, chunk_size, output_dir)

                    csv_file = os.path.join(output_dir, f"genomic_islands_coordinates_{chunk_size}.csv")
                    save_islands_to_fasta_from_csv(genome_seq, csv_file, sequence_id, output_dir, chunk_size // 1000)

                    for info in final_island_info:
                        island_size = info['end_bp'] - info['start_bp']
                        if info['type'] == 'High-Confidence':
                            summary_data["Number of High-Confidence Islands"] += 1
                            summary_data["Total High-Confidence Island Size"] += island_size
                        elif info['type'] == 'Primary':
                            summary_data["Number of Primary Islands"] += 1
                            summary_data["Total Primary Island Size"] += island_size
                        elif info['type'] == 'Secondary':
                            summary_data["Number of Secondary Islands"] += 1
                            summary_data["Total Secondary Island Size"] += island_size

                pd.DataFrame([summary_data]).to_csv(summary_file, mode='a', header=False, index=False)
                print(f"Updated summary for {sequence_id}")

            except Exception as e:
                print(f"FATAL ERROR processing {fasta_file} record {record.id}: {e}")
                import traceback
                traceback.print_exc()
                continue

        if not records_found:
            print(f"Warning: No FASTA records found in {file_path}")

def get_fasta_files_from_folder(folder_path):
    fasta_extensions = (".fasta", ".fa", ".fna")
    return [
        os.path.join(folder_path, name)
        for name in sorted(os.listdir(folder_path))
        if name.lower().endswith(fasta_extensions)
    ]

def process_fasta_folder(folder_path, base_output_dir, min_length, summary_file="genome_summary.csv"):
    fasta_paths = get_fasta_files_from_folder(folder_path)
    if not fasta_paths:
        print(f"Error: No FASTA files found in '{folder_path}'")
        sys.exit(1)
    process_fasta_files(fasta_paths, base_output_dir, min_length, summary_file)

# Helper function to create output directories
def sanitize_sequence_id(sequence_id):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in sequence_id).strip("_")

def create_output_dir(base_output_dir, sequence_id):
    """Creates an output directory for each genome."""
    output_dir = os.path.join(base_output_dir, sequence_id)
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

# Main function to handle command-line arguments
def main():
    """Main entry point of the script."""
    if len(sys.argv) != 3:
        print("Usage: python IslaMine.py <input_fasta_or_folder> <base_output_directory>")
        sys.exit(1)

    input_path = sys.argv[1]
    base_output_dir = sys.argv[2]
    min_length = 50000  # Set a reasonable minimum genome length

    if not os.path.exists(input_path):
        print(f"Error: Input path not found at '{input_path}'")
        sys.exit(1)

    os.makedirs(base_output_dir, exist_ok=True)

    if os.path.isfile(input_path):
        process_fasta_files([input_path], base_output_dir, min_length)
    elif os.path.isdir(input_path):
        process_fasta_folder(input_path, base_output_dir, min_length)
    else:
        print(f"Error: Input path is not a regular FASTA file or folder: '{input_path}'")
        sys.exit(1)

    print("\n--- Analysis Complete ---")

if __name__ == "__main__":
    main()

