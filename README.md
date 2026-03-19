# ASV-to-MAG Matching Pipeline

This Snakemake pipeline matches 16S ASVs (from DADA2) to Metagenome-Assembled Genomes (MAGs) using multiple approaches to achieve both species-level and strain-level resolution.

## Overview

The pipeline performs:
1. **16S rRNA extraction** from MAGs using Barrnap
2. **BLAST-based matching** of ASVs to MAG 16S sequences (species-level)
3. **Whole-genome BLAST** of ASVs to MAG contigs (strain-level)
4. **Marker gene extraction** (optional, for protein-level matching)
5. **Assignment processing** with configurable thresholds
6. **Target species filtering** for specific taxa of interest
7. **Summary statistics** and visualization reports

## Input Format

### ASV Table (`data/asv_table.tsv`)
```tsv
asv_id	asv_seq
ASV_001	/path/to/ASV_001.fasta
ASV_002	/path/to/ASV_002.fasta
```

### MAG Table (`data/mag_table.tsv`)
```tsv
mag_id	mag_seq
MAG_001	/path/to/MAG_001.fasta
MAG_002	/path/to/MAG_002.fasta
```

## Installation

### Required Tools
- Python 3.8+
- Snakemake 6.0+
- BLAST+ (makeblastdb, blastn)
- Barrnap (for 16S extraction)
- DIAMOND (optional, for protein-level matching)

### Install Dependencies
```bash
# Create conda environment
conda create -n asv_to_mag python=3.9 snakemake blast barrnap diamond
conda activate asv_to_mag

# Install Python dependencies
pip install pandas numpy
```

## Configuration

Edit `config.yaml` to customize:

```yaml
# Input files
asv_table: "data/asv_table.tsv"
mag_table: "data/mag_table.tsv"

# Column names for input tables
asv_id_column: "asv_id"
asv_seq_column: "asv_seq"
mag_id_column: "mag_id"
mag_seq_column: "mag_seq"

# Thresholds
identity_threshold: 97  # % identity for species-level
coverage_threshold: 80  # % coverage
strain_threshold: 98.5  # % identity for strain-level

# Target species for filtering
target_species:
  - "Bifidobacterium longum"
  - "Gemella taiwanensis"
  - "Streptococcus oralis"
  - "Klebsiella pneumoniae"

# Computational resources
threads: 16
memory: "64G"
```

## Usage

### Run the Full Pipeline
```bash
snakemake --cores 16 --configfile config.yaml
```

### Run Specific Rules
```bash
# Just extract 16S sequences
snakemake extract_16s_from_mags --cores 16

# Run BLAST only
snakemake blast_asv_to_mag_16s --cores 16

# Generate summary only
snakemake generate_summary --cores 16
```

### Dry Run
```bash
snakemake --cores 16 --dry-run  # See what would be executed
```

## Output Files

```
assignments/
├── asv_to_mag_assignments.tsv    # All ASV-to-MAG assignments
└── target_species_assignments.tsv # Filtered for target species

summary/
├── asv_mag_summary.tsv           # Summary statistics
└── visualization_report.html     # HTML report

blast_results/
├── asv_to_mag_16s.tsv            # 16S BLAST results
└── asv_to_mag_contigs.tsv        # Contig BLAST results
```

## Understanding Results

### Assignment File (`assignments/asv_to_mag_assignments.tsv`)
| Column | Description |
|--------|-------------|
| asv_id | ASV identifier |
| mag_id | MAG identifier |
| identity | Sequence identity (%) |
| alignment_length | Alignment length (bp) |
| evalue | BLAST e-value |
| source | Source of match (16s/contig) |

### Strain Resolution (`assignments/asv_strain_resolution.tsv`)
MAGs with multiple ASV assignments indicate strain-level variation:
- Multiple ASVs per MAG = potential strain variants
- High identity (>98.5%) = same species, different strains

## Workflow Diagram

```
Input Tables
    │
    ├─► Validate Inputs
    │
    ├─► Extract 16S from MAGs (Barrnap)
    │       │
    │       └─► Build 16S BLAST DB
    │               │
    │               └─► BLAST ASVs → 16S Results
    │
    ├─► Build MAG Contig DB
    │       │
    │       └─► BLAST ASVs → Contig Results
    │
    └─► Process BLAST Results
            │
            ├─► Filter by Thresholds
            ├─► Assign Best Hits
            └─► Identify Strain Variants
                    │
                    ├─► Filter Target Species
                    └─► Generate Summary Report
```

## Troubleshooting

### Common Issues

1. **Barrnap not found**
   ```bash
   conda install -c bioconda barrnap
   ```

2. **BLAST database build fails**
   - Check that MAG FASTA files are valid
   - Ensure sufficient disk space

3. **No assignments found**
   - Lower identity threshold in `config.yaml`
   - Verify ASV sequences are full-length 16S
   - Check that MAGs contain 16S regions

## Citation

For questions about the pipeline, please refer to the original tools:
- DADA2: Callahan et al. (2016)
- Barrnap: https://github.com/tseemann/barrnap
- BLAST+: Camacho et al. (2009)

## License

MIT License
