# ATMPI: ASV-to-MAG Integration Pipeline

## 1. Biological & Computational Problem
**The Problem:**
In microbiome research, we often have two powerful but disconnected datasets:
*   **16S rRNA Amplicon (ASVs):** High-resolution taxonomic fingerprints that tell us *who* is there, but not what they are doing.
*   **Metagenomics (MAGs):** Draft genomes that tell us the functional potential (*what* they can do), but often lack 16S markers because these genes are highly repetitive and difficult to assemble or bin correctly.

**The Solution:**
ATMPI bridges this gap by systematically matching ASVs to MAGs. It solves the "missing 16S" problem in metagenomics by using a dual-search strategy: finding matches in predicted 16S genes AND searching the entire genomic content of the MAGs.

## 2. Quick Start: How to Run

### **Step 1: Environment Setup**
Ensure you have Conda or Mamba installed. Activate the project-specific environment:
```bash
# Recommended: Using mamba for speed
mamba activate env-atmpi

# Alternatively: Using conda
conda activate env-atmpi
```

### **Step 2: Prepare Input Tables**
The pipeline requires two tab-separated files: `asv_table.tsv` and `mag_table.tsv`.

**asv_table.tsv format:**
| asv_id | asv_seq |
| :--- | :--- |
| ASV_1 | /path/to/ASV_1.fasta |
| ASV_2 | /path/to/ASV_2.fasta |

**mag_table.tsv format:**
| mag_id | mag_seq |
| :--- | :--- |
| MAG_A | /path/to/MAG_A.fa.gz |
| MAG_B | /path/to/MAG_B.fa.gz |

*Note: MAG files can be gzipped (`.gz`) or plain text.*

### **Step 3: Configuration**
Review `config.yaml` to ensure the column names match your TSV files and adjust thresholds if necessary.

### **Step 4: Execute the Pipeline**
Run the pipeline using Snakemake. It is recommended to perform a dry-run first.

```bash
# 1. Dry-run (checks logic and file paths without running)
snakemake -n --cores 16 --configfile config.yaml

# 2. Real execution (using 16 cores)
snakemake --cores 16 --configfile config.yaml

# 3. Force rerun of a specific part (e.g., if you changed filtering logic)
snakemake --cores 16 --configfile config.yaml -R process_blast_results
```

## 3. Technical Workflow (Rule-by-Rule)

| Step | Rule | Tool | Action |
| :--- | :--- | :--- | :--- |
| **0** | `validate_inputs` | `Python` | Ensures your ASV/MAG tables are formatted correctly and all sequence files exist. |
| **1** | `extract_16s_from_mags` | `Barrnap` | Scans MAGs using HMMs to find ribosomal RNA genes (5S, 16S, 23S). |
| **2** | `build_mag_16s_db` | `makeblastdb` | Indexes the predicted 16S genes for high-speed searching. |
| **3** | `blast_asv_to_mag_16s` | `blastn` | Matches ASVs against the predicted rRNA genes (The "Gold Standard" match). |
| **4** | `blast_asv_to_mag_contigs` | `blastn` | Matches ASVs against **all contigs** in the MAG. This finds 16S sequences that were assembled but not "recognized" by gene predictors. |
| **5** | `process_blast_results` | `Pandas` | **The Brain:** Merges both search results, prioritizes the best hits, and applies filters (97% identity, 80% ASV coverage). |
| **6** | `filter_target_species` | `Python` | Subsets results for specific taxa of interest defined in your config. |
| **7** | `generate_summary` | `Python` | Produces an HTML dashboard and TSV statistics. |

## 3. Understanding the Outputs

### `assignments/asv_to_mag_assignments.tsv`
This is the most important file. Each row represents a high-confidence link between an ASV and a MAG.
*   **identity:** Similarity percentage. >97% is typical for species, >99% for strains.
*   **source:** 
    *   `16s`: The ASV matched a predicted rRNA gene.
    *   `contig`: The ASV matched a raw genomic fragment (essential for fragmented MAGs).

### `assignments/asv_strain_resolution.tsv`
Identifies MAGs that matched **multiple different ASVs**. 
*   **Meaning:** This suggests the MAG contains multiple 16S operons or represents a population of closely related strains (microdiversity) that were binned together.

## 4. Analysis & Visualization in R

Once the pipeline finished, you can load the results into R to link taxonomy to function.

### Loading and Plotting Identity
```R
library(ggplot2)
library(readr)

# Load assignments
mapping <- read_tsv("assignments/asv_to_mag_assignments.tsv")

# Plot the distribution of identity scores
ggplot(mapping, aes(x=identity, fill=source)) +
  geom_histogram(binwidth=0.2, alpha=0.8) +
  theme_minimal() +
  labs(title="ASV-MAG Match Quality", x="Identity (%)", y="Count")
```

### Linking to Abundance
If you have an ASV abundance table (`asv_counts.tsv`), you can see how much of your community is captured in your MAGs:
```R
counts <- read_tsv("asv_counts.tsv")
merged <- inner_join(mapping, counts, by="asv_id")

# Summarize abundance by MAG
mag_abundance <- merged %>% 
  group_by(mag_id) %>% 
  summarise(Total_RelAbund = sum(RelativeAbundance))
```

## 5. Finding Matches for Specific ASVs
If you are interested in a specific ASV (e.g., a biomarker), you can query the results directly:

**Via Terminal:**
```bash
# Search for ASV_1003 and sort by best identity
grep "ASV_1003" assignments/asv_to_mag_assignments.tsv | sort -k3,3rn
```

**Understanding the match:**
1.  Check the **identity**: Is it high enough to trust the assignment?
2.  Check the **source**: If it's `16s`, the match is very robust.
3.  Check **coverage**: In `blast_results/`, you can see if the entire ASV length aligned.

## 6. Configuration
Edit `config.yaml` to adjust:
*   `identity_threshold`: Default 97. Increase to 99 for strain-level matching.
*   `target_species`: List species you want to specifically track.
*   `threads`: Number of CPU cores to use.
