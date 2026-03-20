#!/usr/bin/env python3
"""
Filter ASV-MAG assignments for target species
"""
import argparse
import pandas as pd
import os


def normalize_species_name(name):
    """Normalize species name for matching"""
    if pd.isna(name):
        return ""
    # Convert to lowercase and remove special characters
    name = str(name).lower().strip()
    name = name.replace('_', ' ').replace('"', '').replace("'", "")
    return name


def match_species(assignment_mag, target_species_list):
    """Check if MAG matches any target species"""
    normalized_assignment = normalize_species_name(assignment_mag)

    for target in target_species_list:
        normalized_target = normalize_species_name(target)
        # Check for partial match
        if normalized_target in normalized_assignment or normalized_assignment in normalized_target:
            return target

    return None


def main():
    parser = argparse.ArgumentParser(description='Filter assignments for target species')
    parser.add_argument('--assignments', required=True, help='ASV-MAG assignments')
    parser.add_argument('--mag-table', required=True, help='MAG table with taxonomy')
    parser.add_argument('--target-species', required=True, help='Target species list (comma separated)')
    parser.add_argument('--output', required=True, help='Output file')

    args = parser.parse_args()

    # Split target species
    target_species_list = [s.strip() for s in args.target_species.split(',')]

    print("=" * 60)
    print("FILTERING FOR TARGET SPECIES")
    print("=" * 60)
    print(f"Target species: {', '.join(target_species_list)}")

    # Load assignments
    assignments = pd.read_csv(args.assignments, sep='\t')
    print(f"Total assignments: {len(assignments)}")

    # Load MAG table to get taxonomy
    mag_df = pd.read_csv(args.mag_table, sep='\t')

    # Merge with MAG taxonomy
    assignments = assignments.merge(
        mag_df,
        left_on='mag_id',
        right_on='mag_id',
        how='left'
    )

    # Find species column (try common names)
    species_col = None
    for col in ['species', 'gtid', 'taxon', 'taxonomy']:
        if col in assignments.columns:
            species_col = col
            break

    if species_col is None:
        print("WARNING: No species column found in MAG table")
        print(f"Available columns: {list(assignments.columns)}")
        # Try to use mag_id as fallback
        assignments['species'] = assignments['mag_id']
        species_col = 'species'

    print(f"Using species column: {species_col}")

    # Filter for target species
    assignments['matched_species'] = assignments[species_col].apply(
        lambda x: match_species(x, target_species_list)
    )

    filtered = assignments[assignments['matched_species'].notna()].copy()

    # Rename matched species column
    filtered = filtered.rename(columns={'matched_species': 'target_species'})

    # Save output
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    filtered.to_csv(args.output, sep='\t', index=False)

    # Summary
    print("=" * 60)
    print("FILTERING RESULTS")
    print("=" * 60)
    print(f"  - Total assignments: {len(assignments)}")
    print(f"  - Target species matches: {len(filtered)}")
    print(f"  - Unique MAGs: {filtered['mag_id'].nunique()}")
    print(f"  - Unique ASVs: {filtered['asv_id'].nunique()}")

    # Per-species breakdown
    print("\nPer-species breakdown:")
    species_counts = filtered['target_species'].value_counts()
    for species, count in species_counts.items():
        print(f"  - {species}: {count} ASVs")

    print(f"\n  - Output: {args.output}")
    print("=" * 60)


if __name__ == "__main__":
    main()
