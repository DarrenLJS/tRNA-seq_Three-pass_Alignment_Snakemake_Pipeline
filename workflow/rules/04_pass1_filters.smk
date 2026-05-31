# =============================================================================
# workflow/rules/04_pass1_filters.smk
# Apply two pysam-based filters to the mim-tRNAseq BAMs to identify reads
# from functional, mature tRNAs (proposal Section 3.4):
#
#   Filter 1 — CCA filter
#       The 3′ end of a mature tRNA carries the non-encoded CCA tail.
#       R2 reads in the paired-end library prime from this 3′ end, so the
#       first 3 bases of the aligned R2 query sequence should be "CCA".
#       Reads passing this filter are from fully matured tRNAs.
#
#   Filter 2 — Anticodon concordance filter
#       For each read, the anticodon inferred from the reference locus it
#       maps to must match the GtRNAdb-annotated anticodon.
#       This removes chimeric or mis-mapped reads.
#
# Outputs per sample:
#   {sample}.functional.bam   — reads passing BOTH filters
#   {sample}.filter_stats.tsv — QC metrics (pass rates, used for QC summary)
#
# The CCA+anticodon pass rate is the PRIMARY internal QC metric
# (target: ≥70% of Pass 1-aligned reads; proposal Section 4).
# =============================================================================

rule cca_anticodon_filter:
    """
    Apply CCA and anticodon concordance filters to a mim-tRNAseq BAM.
    Calls workflow/scripts/cca_anticodon_filter.py.
    """
    input:
        bam          = f"{SCRATCH}/pass1_mimtrnaseq/{{sample}}/{{sample}}.bam",
        bai          = f"{SCRATCH}/pass1_mimtrnaseq/{{sample}}/{{sample}}.bam.bai",
        anticodon_map = config["references"]["anticodon_map"],
    output:
        filtered_bam = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.functional.bam",
        filtered_bai = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.functional.bam.bai",
        stats        = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.filter_stats.tsv",
        unmapped_r1  = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.unmapped_R1.fq.gz",
        unmapped_r2  = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.unmapped_R2.fq.gz",
    params:
        outdir       = f"{SCRATCH}/pass1_filters/{{sample}}",
        min_cca_qual = 20,   # min base quality at CCA positions to trust the call
    log:
        f"{SCRATCH}/logs/04_pass1_filters/{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/04_pass1_filters/{{sample}}.tsv",
    conda:
        "../../envs/environment.yaml"
    script:
        "../scripts/cca_anticodon_filter.py"
