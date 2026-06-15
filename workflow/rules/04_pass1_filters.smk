# =============================================================================
# workflow/rules/04_pass1_filters.smk
# Apply two pysam-based filters to the mim-tRNAseq BAMs to identify reads
# from functional, mature tRNAs (proposal Section 3.4):
#
#   Filter 1 — CCA filter
#       The 3′ end of a mature tRNA carries the non-encoded CCA tail.
#       R2 reads prime from this 3′ end in REVERSE COMPLEMENT orientation,
#       so the first 3 bases of R2 should be TGG (= RC of CCA), not CCA.
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
#
# NOTE: mim-tRNAseq aligns R1 only and produces single-end (unpaired) BAMs
# with 100% mapping rate — unmapped reads are never written to the BAM.
# FIX (2026-06-14): trimmed_r1 is now also required so the filter script can
# recover unmapped R1 reads by streaming the FASTQ and writing reads whose
# names are absent from the BAM's mapped set (Bug 2 fix). trimmed_r2 is
# required for the TGG/NGG CCA check (RC orientation) on mapped reads (Bug 1).
# =============================================================================

rule cca_anticodon_filter:
    """
    Apply CCA and anticodon concordance filters to a mim-tRNAseq BAM.
    Calls workflow/scripts/cca_anticodon_filter.py.

    mim-tRNAseq produces single-end BAMs (R1 only, 100% mapping rate).
    trimmed_r1 recovers unmapped reads absent from the BAM (Bug 2 fix).
    trimmed_r2 is used for TGG/NGG CCA check in RC orientation (Bug 1 fix).
    """
    input:
        bam           = f"{SCRATCH}/pass1_mimtrnaseq/{{sample}}/{{sample}}.bam",
        bai           = f"{SCRATCH}/pass1_mimtrnaseq/{{sample}}/{{sample}}.bam.bai",
        anticodon_map = config["references"]["anticodon_map"],
        trimmed_r1    = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_val_1.fq.gz",  # FIX (Bug 2): recover unmapped R1 reads absent from BAM
        trimmed_r2    = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_val_2.fq.gz",  # FIX (Bug 1): CCA check in RC orientation (TGG not CCA)
    output:
        filtered_bam = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.functional.bam",
        filtered_bai = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.functional.bam.bai",
        stats        = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.filter_stats.tsv",
        unmapped_r1  = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.unmapped_R1.fq.gz",
        unmapped_r2  = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.unmapped_R2.fq.gz",
    params:
        outdir       = f"{SCRATCH}/pass1_filters/{{sample}}",
        min_cca_qual = 20,   # min base quality at TGG positions 1+2 (pos 0 skipped — often N on R2)
    log:
        f"{SCRATCH}/logs/04_pass1_filters/{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/04_pass1_filters/{{sample}}.tsv",
    threads: 1
    resources:
        runtime   = 240,
        sge_extra = "-V -l h_vmem=64000M"
    conda:
        "../../envs/environment.yaml"
    script:
        "../scripts/cca_anticodon_filter.py"
