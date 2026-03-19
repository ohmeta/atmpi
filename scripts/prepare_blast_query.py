#!/usr/bin/env python3
"""
Prepare ASV sequences for BLAST from table file
"""
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Prepare ASV sequences for BLAST')
    parser.add_argument('--asv-table', required=True, help='ASV table path')
    parser.add_argument('--seq-column', default='asv_seq', help='Sequence column')
    parser.add_argument('--asv-id-column', default='asv_id', help='ASV ID column')
    parser.add_argument('--output', required=True, help='Output FASTA file')

    args = parser.parse_args()

    print("=" * 60)
    print("PREPARING ASV SEQUENCES FOR BLAST")
    print("=" * 60)

    # Read ASV table
    df = pd.read_csv(args.asv_table, sep='\t')

    # Create output directory
    output_dir = '/'.join(args.output.split('/')[:-1])
    if output_dir:
        import os
        os.makedirs(output_dir, exist_ok=True)

    # Write FASTA file
    with open(args.output, 'w') as f:
        for idx, row in df.iterrows():
            asv_id = row[args.asv_id_column]
            seq_path = row[args.seq_column]

            if not seq_path:
                print(f"  WARNING: Empty sequence for {asv_id}")
                continue

            # Read sequence from file
            with open(seq_path, 'r') as seq_file:
                seq = seq_file.read().strip().replace('\n', '')

            # Write to output FASTA
            f.write(f">{asv_id}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    print(f"  - Processed {len(df)} ASVs")
    print(f"  - Output: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
