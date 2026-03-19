#!/usr/bin/env python3
"""
Process BLAST results and assign ASVs to MAGs
"""
import argparse
import pandas as pd
import numpy as np
import os
import sys


def load_blast_results(filepath):
    """Load BLAST results from TSV"""
    columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
        'qlen', 'slen'
    ]
    df = pd.read_csv(filepath, sep='\t', header=None, names=columns)
    return df


def calculate_coverage(df, asv_lengths, mag_lengths):
    """Calculate ASV and MAG coverage"""
    df['asv_coverage'] = (df['length'] / asv_lengths[df['qseqid']].values * 100).round(2)
    df['mag_coverage'] = (df['length'] / mag_lengths[df['sseqid']].values * 100).round(2)
    return df


def filter_assignments(df, identity_thresh, coverage_thresh):
    """Filter assignments by identity and coverage thresholds"""
    # Calculate sequence lengths from ASV and MAG tables
    df['asv_coverage'] = (df['length'] / df['qlen'] * 100).round(2)
    df['mag_coverage'] = (df['length'] / df['slen'] * 100).round(2)

    # Filter by thresholds
    filtered = df[
        (df['pident'] >= identity_thresh) &
        (df['asv_coverage'] >= coverage_thresh) &
        (df['mag_coverage'] >= coverage_thresh)
    ].copy()

    return filtered


def assign_best_hit(df):
    """Assign best hit per ASV"""
    # Sort by identity (desc) and evalue (asc)
    df_sorted = df.sort_values(['qseqid', 'pident', 'evalue'],
                                ascending=[True, False, True])

    # Keep best hit per ASV
    best_hits = df_sorted.drop_duplicates(subset='qseqid', keep='first')

    return best_hits


def parse_sseqid(sseqid):
    """Parse MAG ID from sseqid (format: MAG_ID|original_header)"""
    if '|' in str(sseqid):
        return sseqid.split('|')[0]
    return sseqid


def identify_strains(assignments):
    """Identify strain-level variation (multiple ASVs per MAG)"""
    # Group by MAG and count ASVs
    mag_asv_counts = assignments.groupby('mag_id')['asv_id'].nunique().reset_index()
    mag_asv_counts.columns = ['mag_id', 'n_asvs']

    # Find MAGs with multiple ASVs (strain variants)
    strain_variants = mag_asv_counts[mag_asv_counts['n_asvs'] > 1]

    return strain_variants


def main():
    parser = argparse.ArgumentParser(description='Process BLAST results')
    parser.add_argument('--asv-table', required=True, help='ASV table path')
    parser.add_argument('--mag-table', required=True, help='MAG table path')
    parser.add_argument('--blast-16s', required=True, help='16S BLAST results')
    parser.add_argument('--blast-contigs', required=True, help='Contig BLAST results')
    parser.add_argument('--output-assignments', required=True, help='Output assignments')
    parser.add_argument('--output-strain', required=True, help='Output strain resolution')
    parser.add_argument('--identity-threshold', type=float, default=97, help='Identity threshold')
    parser.add_argument('--coverage-threshold', type=float, default=80, help='Coverage threshold')
    parser.add_argument('--strain-threshold', type=float, default=98.5, help='Strain threshold')

    args = parser.parse_args()

    print("=" * 60)
    print("PROCESSING BLAST RESULTS")
    print("=" * 60)

    # Load ASV and MAG tables
    asv_df = pd.read_csv(args.asv_table, sep='\t')
    mag_df = pd.read_csv(args.mag_table, sep='\t')

    # Load BLAST results
    print("Loading 16S BLAST results...")
    blast_16s = load_blast_results(args.blast_16s)
    print(f"  - Total hits: {len(blast_16s)}")

    print("Loading contig BLAST results...")
    blast_contigs = load_blast_results(args.blast_contigs)
    print(f"  - Total hits: {len(blast_contigs)}")

    # Combine BLAST results (prioritize contig hits)
    print("Combining BLAST results...")
    blast_contigs['source'] = 'contig'
    blast_16s['source'] = '16s'

    # Merge: prefer contig hits over 16S hits
    combined = blast_contigs.merge(
        blast_16s[['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'source']],
        on=['qseqid', 'sseqid'],
        how='outer',
        suffixes=('_contig', '_16s')
    )

    # Keep best hit per ASV-MAG pair
    combined['pident'] = combined['pident_contig'].fillna(combined['pident_16s'])
    combined['length'] = combined['length_contig'].fillna(combined['length_16s'])
    combined['evalue'] = combined['evalue_contig'].fillna(combined['evalue_16s'])
    combined['source'] = combined['source_contig'].fillna(combined['source_16s'])

    combined = combined.dropna(subset=['qseqid', 'sseqid'])

    # Filter by thresholds
    print(f"Filtering by identity >= {args.identity_threshold}% and coverage >= {args.coverage_threshold}%...")
    filtered = filter_assignments(
        combined,
        args.identity_threshold,
        args.coverage_threshold
    )
    print(f"  - Assignments after filtering: {len(filtered)}")

    # Extract best hit per ASV
    print("Extracting best hit per ASV...")
    assignments = assign_best_hit(filtered)

    # Parse MAG IDs from sseqid
    assignments['mag_id'] = assignments['sseqid'].apply(parse_sseqid)

    # Rename columns for output
    assignments_out = assignments[['qseqid', 'mag_id', 'pident', 'length', 'evalue', 'source']].copy()
    assignments_out.columns = ['asv_id', 'mag_id', 'identity', 'alignment_length', 'evalue', 'source']

    # Add taxonomic information from MAG table
    mag_taxonomy = mag_df[['mag_id']].copy()
    assignments_out = assignments_out.merge(mag_taxonomy, on='mag_id', how='left')

    # Save assignments
    os.makedirs(os.path.dirname(args.output_assignments), exist_ok=True)
    assignments_out.to_csv(args.output_assignments, sep='\t', index=False)
    print(f"  - Saved assignments: {args.output_assignments}")

    # Identify strain-level variation
    print("Identifying strain-level variation...")
    strain_variants = identify_strains(assignments_out)

    # Save strain resolution
    os.makedirs(os.path.dirname(args.output_strain), exist_ok=True)
    strain_variants.to_csv(args.output_strain, sep='\t', index=False)
    print(f"  - Saved strain resolution: {args.output_strain}")
    print(f"  - MAGs with strain variants: {len(strain_variants)}")

    # Summary statistics
    print("=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    print(f"  - Total ASVs: {asv_df.shape[0]}")
    print(f"  - Total MAGs: {mag_df.shape[0]}")
    print(f"  - ASVs with MAG assignments: {assignments_out['asv_id'].nunique()}")
    print(f"  - MAGs with ASV assignments: {assignments_out['mag_id'].nunique()}")
    print(f"  - MAGs with strain variants: {len(strain_variants)}")
    print("=" * 60)


if __name__ == "__main__":
    main()
