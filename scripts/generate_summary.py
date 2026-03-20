#!/usr/bin/env python3
"""
Generate summary statistics and visualization report
"""
import argparse
import pandas as pd
import numpy as np
import os
from datetime import datetime


def generate_summary_table(assignments_df, output_path):
    """Generate summary statistics table"""
    summary = {
        'metric': [],
        'value': [],
        'description': []
    }

    # Overall statistics
    summary['metric'].append('total_asvs')
    summary['value'].append(assignments_df['asv_id'].nunique())
    summary['description'].append('Total unique ASVs')

    summary['metric'].append('total_mags')
    summary['value'].append(assignments_df['mag_id'].nunique())
    summary['description'].append('Total unique MAGs')

    summary['metric'].append('total_assignments')
    summary['value'].append(len(assignments_df))
    summary['description'].append('Total ASV-MAG assignments')

    # Identity statistics
    summary['metric'].append('mean_identity')
    summary['value'].append(round(assignments_df['identity'].mean(), 2))
    summary['description'].append('Mean identity (%)')

    summary['metric'].append('min_identity')
    summary['value'].append(round(assignments_df['identity'].min(), 2))
    summary['description'].append('Minimum identity (%)')

    summary['metric'].append('max_identity')
    summary['value'].append(round(assignments_df['identity'].max(), 2))
    summary['description'].append('Maximum identity (%)')

    # Source statistics
    source_counts = assignments_df['source'].value_counts()
    for source, count in source_counts.items():
        summary['metric'].append(f'{source}_assignments')
        summary['value'].append(count)
        summary['description'].append(f'{source.capitalize()} assignments')

    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(output_path, sep='\t', index=False)

    return summary_df


def generate_per_species_summary(assignments_df, output_path):
    """Generate per-species breakdown"""
    # Group by MAG and ASV
    per_mag = assignments_df.groupby('mag_id').agg({
        'asv_id': 'nunique',
        'identity': ['mean', 'min', 'max']
    }).reset_index()

    per_mag.columns = ['mag_id', 'n_asvs', 'mean_identity', 'min_identity', 'max_identity']

    # Sort by number of ASVs
    per_mag = per_mag.sort_values('n_asvs', ascending=False)

    per_mag.to_csv(output_path, sep='\t', index=False)
    return per_mag


def generate_html_report(assignments_df, summary_df, output_path):
    """Generate HTML visualization report"""

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>ASV-to-MAG Assignment Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #333;
            border-bottom: 2px solid #007bff;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #555;
            margin-top: 30px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #007bff;
            color: white;
            font-weight: bold;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .metric-card {{
            display: inline-block;
            margin: 10px;
            padding: 20px;
            background-color: #007bff;
            color: white;
            border-radius: 8px;
            text-align: center;
            min-width: 150px;
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
        }}
        .metric-label {{
            font-size: 0.9em;
            opacity: 0.9;
        }}
        .timestamp {{
            color: #666;
            font-size: 0.9em;
            margin-top: 20px;
        }}
        .warning {{
            background-color: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 10px;
            margin: 10px 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>ASV-to-MAG Assignment Report</h1>

        <h2>Summary Metrics</h2>
        <div>
"""

    # Add metric cards
    for _, row in summary_df.iterrows():
        if row['metric'] in ['total_asvs', 'total_mags', 'total_assignments', 'mean_identity']:
            html += f"""
            <div class="metric-card">
                <div class="metric-value">{row['value']}</div>
                <div class="metric-label">{row['description']}</div>
            </div>
"""

    html += """
        </div>

        <h2>Assignment Sources</h2>
        <table>
            <tr>
                <th>Source</th>
                <th>Number of Assignments</th>
                <th>Percentage</th>
            </tr>
"""

    source_counts = assignments_df['source'].value_counts()
    total = len(assignments_df)
    for source, count in source_counts.items():
        pct = (count / total * 100)
        html += f"""
            <tr>
                <td>{source}</td>
                <td>{count}</td>
                <td>{pct:.1f}%</td>
            </tr>
"""

    html += """
        </table>

        <h2>Top MAGs by ASV Assignments</h2>
        <table>
            <tr>
                <th>Rank</th>
                <th>MAG ID</th>
                <th>Number of ASVs</th>
                <th>Avg Identity (%)</th>
                <th>Identity Range (%)</th>
            </tr>
"""

    # Get top 20 MAGs
    top_mags = assignments_df.groupby('mag_id').agg({
        'asv_id': 'nunique',
        'identity': ['mean', 'min', 'max']
    }).reset_index()
    top_mags.columns = ['mag_id', 'n_asvs', 'mean_identity', 'min_identity', 'max_identity']
    top_mags = top_mags.sort_values('n_asvs', ascending=False).head(20)

    for idx, row in top_mags.iterrows():
        html += f"""
            <tr>
                <td>{idx + 1}</td>
                <td>{row['mag_id']}</td>
                <td>{row['n_asvs']}</td>
                <td>{row['mean_identity']:.1f}</td>
                <td>{row['min_identity']:.1f}-{row['max_identity']:.1f}</td>
            </tr>
"""

    html += f"""
        </table>

        <div class="timestamp">
            Report generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </div>
</body>
</html>
"""

    with open(output_path, 'w') as f:
        f.write(html)


def main():
    parser = argparse.ArgumentParser(description='Generate summary statistics')
    parser.add_argument('--assignments', required=True, help='ASV-MAG assignments')
    parser.add_argument('--taxonomy', required=False, help='MAG taxonomy file (optional)')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')

    args = parser.parse_args()

    print("=" * 60)
    print("GENERATING SUMMARY")
    print("=" * 60)

    # Load assignments
    assignments_df = pd.read_csv(args.assignments, sep='\t')
    print(f"Loaded {len(assignments_df)} assignments")

    # Create output directory
    output_dir = args.output_prefix
    if not output_dir.endswith('.tsv'):
        os.makedirs(output_dir, exist_ok=True)

    # Generate summary table
    summary_path = f"{output_dir}/asv_mag_summary.tsv"
    summary_df = generate_summary_table(assignments_df, summary_path)
    print(f"Summary saved: {summary_path}")

    # Generate per-species summary
    per_species_path = f"{output_dir}/per_mag_summary.tsv"
    generate_per_species_summary(assignments_df, per_species_path)
    print(f"Per-MAG summary saved: {per_species_path}")

    # Generate HTML report
    html_path = f"{output_dir}/visualization_report.html"
    generate_html_report(assignments_df, summary_df, html_path)
    print(f"HTML report saved: {html_path}")

    print("=" * 60)
    print("SUMMARY COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
