#!/usr/bin/env python3
"""
Validate input ASV and MAG table files
"""
import argparse
import pandas as pd
import os
import sys


def validate_table(table_path, id_col, seq_col, file_type):
    """Validate a table file and check if sequence files exist"""
    print(f"Validating {file_type} table: {table_path}")

    if not os.path.exists(table_path):
        print(f"ERROR: {file_type} table not found: {table_path}")
        return False

    try:
        df = pd.read_csv(table_path, sep='\t')
    except Exception as e:
        print(f"ERROR: Failed to read {file_type} table: {e}")
        return False

    # Check required columns
    if id_col not in df.columns:
        print(f"ERROR: {id_col} column not found in {file_type} table")
        return False

    if seq_col not in df.columns:
        print(f"ERROR: {seq_col} column not found in {file_type} table")
        return False

    # Count records
    n_records = len(df)
    n_valid_seqs = df[seq_col].apply(lambda x: os.path.exists(x)).sum()

    print(f"  - Records: {n_records}")
    print(f"  - Valid sequence files: {n_valid_seqs}/{n_records}")

    if n_valid_seqs < n_records:
        missing = df[~df[seq_col].apply(os.path.exists)][seq_col].tolist()
        print(f"  - WARNING: {len(missing)} missing sequence files:")
        for m in missing[:5]:  # Show first 5
            print(f"      {m}")

    # Show sample
    print(f"  - Sample {id_col}: {df[id_col].iloc[0] if n_records > 0 else 'N/A'}")
    print(f"  - Sample {seq_col}: {df[seq_col].iloc[0] if n_records > 0 else 'N/A'}")

    return True


def main():
    parser = argparse.ArgumentParser(description='Validate input files')
    parser.add_argument('--asv-table', required=True, help='Path to ASV table')
    parser.add_argument('--mag-table', required=True, help='Path to MAG table')
    parser.add_argument('--asv-id-column', default='asv_id', help='ASV ID column name')
    parser.add_argument('--asv-seq-column', default='asv_seq', help='ASV sequence column name')
    parser.add_argument('--mag-id-column', default='mag_id', help='MAG ID column name')
    parser.add_argument('--mag-seq-column', default='mag_seq', help='MAG sequence column name')

    args = parser.parse_args()

    print("=" * 60)
    print("INPUT VALIDATION")
    print("=" * 60)

    asv_valid = validate_table(
        args.asv_table,
        args.asv_id_column,
        args.asv_seq_column,
        "ASV"
    )

    mag_valid = validate_table(
        args.mag_table,
        args.mag_id_column,
        args.mag_seq_column,
        "MAG"
    )

    print("=" * 60)
    if asv_valid and mag_valid:
        print("VALIDATION PASSED")
        print("=" * 60)
        sys.exit(0)
    else:
        print("VALIDATION FAILED")
        print("=" * 60)
        sys.exit(1)


if __name__ == "__main__":
    main()
