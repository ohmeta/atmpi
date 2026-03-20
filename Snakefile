# Snakemake Pipeline: ASV-to-MAG Matching
# Matches 16S ASVs to MAGs using multiple approaches

configfile: "config.yaml"

# ============================================================================
# MAIN TARGET
# ============================================================================
rule all:
    input:
        "assignments/asv_to_mag_assignments.tsv",
        "assignments/target_species_assignments.tsv",
        "summary/asv_mag_summary.tsv",
        "summary/visualization_report.html"

# ============================================================================
# CONSTANTS
# ============================================================================
BLAST_EVALUE = 0.0001
BLAST_IDENTITY_THRESHOLD = 97  # % for species-level
BLAST_COVERAGE_THRESHOLD = 80  # % for ASV coverage
STRAIN_IDENTITY_THRESHOLD = 98.5  # % for strain-level

# ============================================================================
# RULES
# ============================================================================

# -----------------------------------------------------------------------------
# Rule 0: Validate input files
# -----------------------------------------------------------------------------
rule validate_inputs:
    input:
        asv_table = config["asv_table"],
        mag_table = config["mag_table"]
    output:
        temp("validation/complete.flag")
    params:
        asv_seq_column = config["asv_seq_column"],
        mag_seq_column = config["mag_seq_column"],
        asv_id_column = config["asv_id_column"],
        mag_id_column = config["mag_id_column"]
    shell:
        """
        echo "Validating input files..."
        python scripts/validate_inputs.py \
            --asv-table {input.asv_table} \
            --mag-table {input.mag_table} \
            --asv-seq-column {params.asv_seq_column} \
            --mag-seq-column {params.mag_seq_column} \
            --asv-id-column {params.asv_id_column} \
            --mag-id-column {params.mag_id_column}
        touch {output[0]}
        """

# -----------------------------------------------------------------------------
# Rule 1: Extract 16S rRNA genes from MAGs using Barrnap
# -----------------------------------------------------------------------------
rule extract_16s_from_mags:
    input:
        mag_table = config["mag_table"]
    output:
        outdir = directory("16s_extraction/mag_16s"),
        fasta = "16s_extraction/mag_16s/all_mag_16s.fasta"
    params:
        mag_seq_column = config["mag_seq_column"],
        mag_id_column = config["mag_id_column"]
    threads:
        config["threads"]
    log:
        "logs/16s_extraction.log"
    shell:
        """
        mkdir -p {output.outdir} logs
        echo "Extracting 16S rRNA genes from MAGs..."

        # Process each MAG file
        python scripts/extract_16s.py \
            --mag-table {input.mag_table} \
            --seq-column {params.mag_seq_column} \
            --mag-id-column {params.mag_id_column} \
            --output-dir {output.outdir} \
            --threads {threads} \
            2> {log}

        echo "16S extraction complete. See {log}"
        """

# -----------------------------------------------------------------------------
# Rule 2: Build BLAST database from MAG 16S sequences
# -----------------------------------------------------------------------------
rule build_mag_16s_db:
    input:
        "16s_extraction/mag_16s/all_mag_16s.fasta"
    output:
        multiext("blast_db/mag_16s_db", ".nhr", ".nin", ".nsq")
    params:
        threads = config["threads"]
    log:
        "logs/build_mag_16s_db.log"
    shell:
        """
        echo "Building BLAST database from MAG 16S sequences..."
        makeblastdb \
            -in {input[0]} \
            -dbtype nucl \
            -out blast_db/mag_16s_db \
            -title "MAG_16S_Database" \
            2> {log}
        """

# -----------------------------------------------------------------------------
# Rule 3: BLAST ASVs against MAG 16S sequences (Species-level)
# -----------------------------------------------------------------------------
rule blast_asv_to_mag_16s:
    input:
        asv_table = config["asv_table"],
        mag_db = multiext("blast_db/mag_16s_db", ".nhr", ".nin", ".nsq")
    output:
        "blast_results/asv_to_mag_16s.tsv"
    params:
        evalue = config["blast_evalue"],
        max_target_seqs = config["blast_max_target_seqs"],
        task = config["blast_task"],
        asv_seq_column = config["asv_seq_column"],
        asv_id_column = config["asv_id_column"],
        db_prefix = "blast_db/mag_16s_db"
    threads:
        config["threads"]
    log:
        "logs/blast_asv_to_mag_16s.log"
    shell:
        """
        echo "BLASTing ASVs against MAG 16S sequences..."

        # Prepare ASV sequences
        python scripts/prepare_blast_query.py \
            --asv-table {input.asv_table} \
            --seq-column {params.asv_seq_column} \
            --asv-id-column {params.asv_id_column} \
            --output blast_query/asv_seqs.fasta

        # Run BLAST
        blastn \
            -query blast_query/asv_seqs.fasta \
            -db {params.db_prefix} \
            -out {output[0]} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
            -evalue {params.evalue} \
            -max_target_seqs {params.max_target_seqs} \
            -task {params.task} \
            -num_threads {threads} \
            2> {log}
        """

# -----------------------------------------------------------------------------
# Rule 8: BLAST ASVs against MAG contigs (Whole genome alignment)
# -----------------------------------------------------------------------------
rule blast_asv_to_mag_contigs:
    input:
        asv_table = config["asv_table"],
        mag_table = config["mag_table"]
    output:
        "blast_results/asv_to_mag_contigs.tsv"
    params:
        evalue = config["blast_evalue"],
        max_target_seqs = config["blast_max_target_seqs"],
        mag_seq_column = config["mag_seq_column"],
        mag_id_column = config["mag_id_column"],
        asv_seq_column = config["asv_seq_column"],
        asv_id_column = config["asv_id_column"],
        task = config["blast_task"]
    threads:
        config["threads"]
    log:
        "logs/blast_asv_to_mag_contigs.log"
    shell:
        """
        echo "BLASTing ASVs against MAG contigs..."

        # Build combined MAG contig database
        python scripts/build_mag_contig_db.py \
            --mag-table {input.mag_table} \
            --seq-column {params.mag_seq_column} \
            --mag-id-column {params.mag_id_column} \
            --output blast_db/mag_contigs_db

        # Prepare ASV sequences (already done in other rules if possible, but keep for independence)
        python scripts/prepare_blast_query.py \
            --asv-table {input.asv_table} \
            --seq-column {params.asv_seq_column} \
            --asv-id-column {params.asv_id_column} \
            --output blast_query/asv_seqs.fasta

        # Run BLAST
        blastn \
            -query blast_query/asv_seqs.fasta \
            -db blast_db/mag_contigs_db \
            -out {output[0]} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
            -evalue {params.evalue} \
            -max_target_seqs {params.max_target_seqs} \
            -task {params.task} \
            -num_threads {threads} \
            2> {log}
        """

# -----------------------------------------------------------------------------
# Rule 9: Process BLAST results and assign ASVs to MAGs
# -----------------------------------------------------------------------------
rule process_blast_results:
    input:
        asv_table = config["asv_table"],
        mag_table = config["mag_table"],
        blast_16s = "blast_results/asv_to_mag_16s.tsv",
        blast_contigs = "blast_results/asv_to_mag_contigs.tsv"
    output:
        "assignments/asv_to_mag_assignments.tsv",
        "assignments/asv_strain_resolution.tsv"
    params:
        identity_threshold = config["identity_threshold"],
        coverage_threshold = config["coverage_threshold"],
        strain_threshold = config["strain_threshold"]
    log:
        "logs/process_blast_results.log"
    shell:
        """
        python scripts/process_blast_results.py \
            --asv-table {input.asv_table} \
            --mag-table {input.mag_table} \
            --blast-16s {input.blast_16s} \
            --blast-contigs {input.blast_contigs} \
            --output-assignments {output[0]} \
            --output-strain {output[1]} \
            --identity-threshold {params.identity_threshold} \
            --coverage-threshold {params.coverage_threshold} \
            --strain-threshold {params.strain_threshold} \
            2> {log}
        """

# -----------------------------------------------------------------------------
# Rule 10: Generate summary statistics
# -----------------------------------------------------------------------------
rule generate_summary:
    input:
        assignments = "assignments/asv_to_mag_assignments.tsv",
        mag_taxonomy = config["mag_taxonomy"] if "mag_taxonomy" in config and config["mag_taxonomy"] else []
    output:
        "summary/asv_mag_summary.tsv",
        "summary/visualization_report.html"
    log:
        "logs/generate_summary.log"
    shell:
        """
        TAX_ARG=""
        if [ ! -z "{input.mag_taxonomy}" ]; then
            TAX_ARG="--taxonomy {input.mag_taxonomy}"
        fi

        python scripts/generate_summary.py \
            --assignments {input.assignments} \
            $TAX_ARG \
            --output-prefix summary \
            2> {log}
        """

# -----------------------------------------------------------------------------
# Rule 11: Filter for target species
# -----------------------------------------------------------------------------
rule filter_target_species:
    input:
        assignments = "assignments/asv_to_mag_assignments.tsv",
        mag_table = config["mag_table"]
    output:
        "assignments/target_species_assignments.tsv"
    params:
        target_species = lambda wildcards: ",".join(config["target_species"])
    log:
        "logs/filter_target_species.log"
    shell:
        """
        python scripts/filter_target_species.py \
            --assignments {input.assignments} \
            --mag-table {input.mag_table} \
            --target-species "{params.target_species}" \
            --output {output[0]} \
            2> {log}
        """
