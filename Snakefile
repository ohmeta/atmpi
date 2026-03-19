# Snakemake Pipeline: ASV-to-MAG Matching
# Matches 16S ASVs to MAGs using multiple approaches

configfile: "config.yaml"

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
    shell:
        """
        echo "Validating input files..."
        python scripts/validate_inputs.py \
            --asv-table {input.asv_table} \
            --mag-table {input.mag_table} \
            --asv-seq-column {config.asv_seq_column} \
            --mag-seq-column {config.mag_seq_column} \
            --asv-id-column {config.asv_id_column} \
            --mag-id-column {config.mag_id_column}
        touch {output[0]}
        """

# -----------------------------------------------------------------------------
# Rule 1: Extract 16S rRNA genes from MAGs using Barrnap
# -----------------------------------------------------------------------------
rule extract_16s_from_mags:
    input:
        mag_table = config["mag_table"]
    output:
        dir("16s_extraction/mag_16s")
    params:
        threads = config["threads"]
    log:
        "logs/16s_extraction.log"
    shell:
        """
        mkdir -p {output.dir} logs
        echo "Extracting 16S rRNA genes from MAGs..."

        # Process each MAG file
        python scripts/extract_16s.py \
            --mag-table {input.mag_table} \
            --seq-column {config.mag_seq_column} \
            --mag-id-column {config.mag_id_column} \
            --output-dir {output.dir} \
            --threads {params.threads}

        echo "16S extraction complete. See {params.log}"
        """

# -----------------------------------------------------------------------------
# Rule 2: Build BLAST database from MAG 16S sequences
# -----------------------------------------------------------------------------
rule build_mag_16s_db:
    input:
        "16s_extraction/mag_16s/all_mag_16s.fasta"
    output:
        "blast_db/mag_16s_db.{n,2.nhr,2.nin,2.nsq}"
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
            -parse_seqids \
            -title "MAG_16S_Database" \
            2> {params.log}
        """

# -----------------------------------------------------------------------------
# Rule 3: BLAST ASVs against MAG 16S sequences (Species-level)
# -----------------------------------------------------------------------------
rule blast_asv_to_mag_16s:
    input:
        asv_table = config["asv_table"],
        mag_db = "blast_db/mag_16s_db"
    output:
        "blast_results/asv_to_mag_16s.tsv"
    params:
        evalue = config["blast_evalue"],
        max_target_seqs = config["blast_max_target_seqs"],
        threads = config["threads"],
        task = "blastn"  # Use "blastn-short" for short reads
    log:
        "logs/blast_asv_to_mag_16s.log"
    shell:
        """
        echo "BLASTing ASVs against MAG 16S sequences..."

        # Prepare ASV sequences
        python scripts/prepare_blast_query.py \
            --asv-table {input.asv_table} \
            --seq-column {config.asv_seq_column} \
            --asv-id-column {config.asv_id_column} \
            --output blast_query/asv_seqs.fasta

        # Run BLAST
        blastn \
            -query blast_query/asv_seqs.fasta \
            -db {input.mag_db} \
            -out {output[0]} \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
            -evalue {params.evalue} \
            -max_target_seqs {params.max_target_seqs} \
            -task {params.task} \
            -num_threads {params.threads} \
            2> {params.log}
        """

# -----------------------------------------------------------------------------
# Rule 4: Extract marker genes from MAGs (GTDB-Tk markers)
# -----------------------------------------------------------------------------
rule extract_marker_genes:
    input:
        mag_table = config["mag_table"]
    output:
        dir("marker_extraction/marker_genes")
    params:
        threads = config["threads"]
    log:
        "logs/marker_extraction.log"
    shell:
        """
        echo "Extracting marker genes from MAGs..."
        mkdir -p {output.dir} logs

        python scripts/extract_markers.py \
            --mag-table {input.mag_table} \
            --seq-column {config.mag_seq_column} \
            --mag-id-column {config.mag_id_column} \
            --output-dir {output.dir} \
            --marker-type {config.marker_type} \
            --threads {params.threads}

        echo "Marker extraction complete."
        """

# -----------------------------------------------------------------------------
# Rule 5: Build DIAMOND database from MAG marker proteins
# -----------------------------------------------------------------------------
rule build_marker_db:
    input:
        "marker_extraction/marker_genes/all_markers.faa"
    output:
        "diamond_db/marker_db.dmnd"
    params:
        threads = config["threads"]
    log:
        "logs/build_marker_db.log"
    shell:
        """
        echo "Building DIAMOND database from marker proteins..."
        diamond makedb \
            --in {input[0]} \
            --db diamond_db/marker_db \
            --threads {params.threads} \
            2> {params.log}
        """

# -----------------------------------------------------------------------------
# Rule 6: Extract protein sequences from ASVs (if needed)
# -----------------------------------------------------------------------------
rule extract_asv_proteins:
    input:
        asv_table = config["asv_table"]
    output:
        "protein_extraction/asv_proteins.faa"
    params:
        frame = 1
    log:
        "logs/extract_asv_proteins.log"
    shell:
        """
        echo "Extracting protein sequences from ASVs..."
        mkdir -p protein_extraction logs

        # Translate ASV nucleotide sequences to proteins
        python scripts/translate_asvs.py \
            --asv-table {input.asv_table} \
            --seq-column {config.asv_seq_column} \
            --asv-id-column {config.asv_id_column} \
            --output {output[0]} \
            --frame {params.frame}

        echo "Protein extraction complete."
        """

# -----------------------------------------------------------------------------
# Rule 7: DIAMOND BLAST ASVs against MAG marker proteins
# -----------------------------------------------------------------------------
rule diamond_asv_to_mag_markers:
    input:
        asv_proteins = "protein_extraction/asv_proteins.faa",
        mag_markers = "diamond_db/marker_db.dmnd"
    output:
        "blast_results/asv_to_mag_markers.tsv"
    params:
        evalue = config["blast_evalue"],
        max_target_seqs = config["blast_max_target_seqs"],
        threads = config["threads"]
    log:
        "logs/diamond_asv_to_mag_markers.log"
    shell:
        """
        echo "Running DIAMOND BLAST on ASVs vs MAG markers..."

        diamond blastp \
            -q {input.asv_proteins} \
            -d {input.mag_markers} \
            -o {output[0]} \
            -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
            -e {params.evalue} \
            -k {params.max_target_seqs} \
            --threads {params.threads} \
            2> {params.log}
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
        threads = config["threads"]
    log:
        "logs/blast_asv_to_mag_contigs.log"
    shell:
        """
        echo "BLASTing ASVs against MAG contigs..."

        # Build combined MAG contig database
        python scripts/build_mag_contig_db.py \
            --mag-table {input.mag_table} \
            --seq-column {config.mag_seq_column} \
            --mag-id-column {config.mag_id_column} \
            --output blast_db/mag_contigs_db

        # Prepare ASV sequences
        python scripts/prepare_blast_query.py \
            --asv-table {input.asv_table} \
            --seq-column {config.asv_seq_column} \
            --asv-id-column {config.asv_id_column} \
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
            -num_threads {params.threads} \
            2> {params.log}
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
    run:
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
                2> {params.log}
            """

# -----------------------------------------------------------------------------
# Rule 10: Generate summary statistics
# -----------------------------------------------------------------------------
rule generate_summary:
    input:
        assignments = "assignments/asv_to_mag_assignments.tsv",
        mag_taxonomy = config["mag_taxonomy"] if "mag_taxonomy" in config else None
    output:
        "summary/asv_mag_summary.tsv",
        "summary/visualization_report.html"
    log:
        "logs/generate_summary.log"
    run:
        shell:
            """
            python scripts/generate_summary.py \
                --assignments {input.assignments} \
                --taxonomy {input.mag_taxonomy} \
                --output-prefix summary \
                2> {params.log}
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
        target_species = config["target_species"]
    log:
        "logs/filter_target_species.log"
    run:
        shell:
            """
            python scripts/filter_target_species.py \
                --assignments {input.assignments} \
                --mag-table {input.mag_table} \
                --target-species {params.target_species} \
                --output {output[0]} \
                2> {params.log}
            """

# ============================================================================
# MAIN TARGET
# ============================================================================
rule all:
    input:
        "assignments/asv_to_mag_assignments.tsv",
        "assignments/target_species_assignments.tsv",
        "summary/asv_mag_summary.tsv",
        "summary/visualization_report.html"
