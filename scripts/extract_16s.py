#!/usr/bin/env python3
"""
Extract 16S rRNA genes from MAGs using Barrnap
"""
import argparse
import pandas as pd
import os
import subprocess
import tempfile
from pathlib import Path


def extract_16s_barrnap(fasta_path, output_dir, threads=8):
    """Extract 16S rRNA genes using Barrnap"""
    print(f"  Extracting 16S from: {fasta_path}")

    # Create temp file for Barrnap output
    with tempfile.NamedTemporaryFile(mode='w', suffix='_16s.fasta',
                                     dir=output_dir, delete=False) as tmp:
        tmp_path = tmp.name

    try:
        cmd = [
            'barrnap',
            '--threads', str(threads),
            '--genomic',
            fasta_path,
            tmp_path
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            print(f"    WARNING: Barrnap failed for {fasta_path}")
            print(f"    STDERR: {result.stderr}")
            return None

        # Read and parse output
        if os.path.exists(tmp_path):
            with open(tmp_path, 'r') as f:
                return f.read()
        return None

    except subprocess.TimeoutExpired:
        print(f"    ERROR: Barrnap timeout for {fasta_path}")
        return None
    except Exception as e:
        print(f"    ERROR: {e}")
        return None
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def read_fasta_sequences(fasta_path):
    """Read a FASTA file and return dict of sequences"""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def write_fasta(sequences, output_path):
    """Write sequences to FASTA file"""
    with open(output_path, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")


def main():
    parser = argparse.ArgumentParser(description='Extract 16S from MAGs')
    parser.add_argument('--mag-table', required=True, help='MAG table path')
    parser.add_argument('--seq-column', default='mag_seq', help='Sequence column')
    parser.add_argument('--mag-id-column', default='mag_id', help='MAG ID column')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=8, help='Number of threads')

    args = parser.parse_args()

    print("=" * 60)
    print("16S EXTRACTION FROM MAGS")
    print("=" * 60)

    # Read MAG table
    df = pd.read_csv(args.mag_table, sep='\t')

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Collect all 16S sequences
    all_16s = {}
    extraction_stats = {"success": 0, "failed": 0}

    for idx, row in df.iterrows():
        mag_id = row[args.mag_id_column]
        seq_path = row[args.seq_column]

        if not os.path.exists(seq_path):
            print(f"  Skipping {mag_id}: file not found")
            extraction_stats["failed"] += 1
            continue

        # Extract 16S
        extracted = extract_16s_barrnap(seq_path, args.output_dir, args.threads)

        if extracted:
            # Parse extracted 16S sequences
            extracted_seqs = read_fasta_sequences(extracted)
            for seq_id, seq in extracted_seqs.items():
                # Prefix with MAG ID
                new_id = f"{mag_id}_{seq_id}"
                all_16s[new_id] = seq
                extraction_stats["success"] += 1
        else:
            extraction_stats["failed"] += 1

    # Write combined 16S file
    output_path = os.path.join(args.output_dir, "all_mag_16s.fasta")
    write_fasta(all_16s, output_path)

    print("=" * 60)
    print(f"EXTRACTION COMPLETE")
    print(f"  - Successful: {extraction_stats['success']}")
    print(f"  - Failed: {extraction_stats['failed']}")
    print(f"  - Total sequences: {len(all_16s)}")
    print(f"  - Output: {output_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
