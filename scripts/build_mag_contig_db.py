#!/usr/bin/env python3
"""
Build combined MAG contig database for BLAST
"""
import argparse
import pandas as pd
import os
import subprocess


def main():
    parser = argparse.ArgumentParser(description='Build MAG contig database for BLAST')
    parser.add_argument('--mag-table', required=True, help='MAG table path')
    parser.add_argument('--seq-column', default='mag_seq', help='Sequence column')
    parser.add_argument('--mag-id-column', default='mag_id', help='MAG ID column')
    parser.add_argument('--output', required=True, help='Output database prefix')

    args = parser.parse_args()

    print("=" * 60)
    print("BUILDING MAG CONTIG DATABASE")
    print("=" * 60)

    # Read MAG table
    df = pd.read_csv(args.mag_table, sep='\t')

    # Create output directory
    output_dir = '/'.join(args.output.split('/')[:-1])
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # Combine all MAG contigs into single file
    combined_path = f"{output_dir}/combined_mag_contigs.fasta"
    import gzip
    with open(combined_path, 'w') as outfile:
        for idx, row in df.iterrows():
            mag_id = row[args.mag_id_column]
            seq_path = row[args.seq_column]

            if not os.path.exists(seq_path):
                print(f"  WARNING: File not found: {seq_path}")
                continue

            # Read and rewrite with modified headers
            if str(seq_path).endswith('.gz'):
                infile = gzip.open(seq_path, 'rt')
            else:
                infile = open(seq_path, 'r')
            
            with infile:
                for line in infile:
                    if line.startswith('>'):
                        # Modify header to include MAG ID
                        header = line[1:].strip()
                        outfile.write(f">{mag_id}|{header}\n")
                    else:
                        outfile.write(line)

    print(f"  - Combined {len(df)} MAG files")
    print(f"  - Combined file: {combined_path}")

    # Build BLAST database
    db_prefix = f"{output_dir}/mag_contigs_db"
    cmd = [
        'makeblastdb',
        '-in', combined_path,
        '-dbtype', 'nucl',
        '-out', db_prefix,
        '-title', 'MAG_Contigs_Database'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: BLAST database build failed")
        print(f"  STDERR: {result.stderr}")
        return

    print(f"  - BLAST database: {db_prefix}")
    print("=" * 60)


if __name__ == "__main__":
    main()
