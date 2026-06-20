# =============================================================================
# workflow/rules/08_deseq2_prep.smk
# Collate mim-tRNAseq count tables into per-cell-line count matrices and
# write the corresponding coldata (sample metadata) TSVs for DESeq2/edgeR.
#
# Also assembles a per-cell-line read assignment summary (QC metric 2 from
# the proposal: proportional read assignment across all categories).
#
# Outputs per cell line:
#   isodecoder_counts_matrix.tsv   — rows=isodecoders, cols=samples
#   isoacceptor_counts_matrix.tsv  — rows=isoacceptors, cols=samples
#   coldata.tsv                    — sample metadata for DESeq2
#   read_assignment_summary.tsv    — read proportions per sample (QC)
# =============================================================================

rule build_deseq2_inputs:
    """
    Collate per-cell-line count matrices and coldata for DESeq2.
    Calls workflow/scripts/collate_counts.py.

    FIX: manifest path is now explicitly passed as params.manifest so that
    collate_counts.py can read sample metadata (condition, replicate, etc.)
    to build coldata.tsv. Access it in the script via snakemake.params.manifest.
    """
    input:
        iso_counts   = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isodecoder_counts.txt",
        isoa_counts  = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isoacceptor_counts.txt",
    output:
        iso_matrix   = f"{SCRATCH}/deseq2_input/{{cell_line}}/isodecoder_counts_matrix.tsv",
        isoa_matrix  = f"{SCRATCH}/deseq2_input/{{cell_line}}/isoacceptor_counts_matrix.tsv",
        coldata      = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    params:
        outdir       = f"{SCRATCH}/deseq2_input/{{cell_line}}",
        cell_line    = "{cell_line}",
        manifest     = config["manifest"],   # FIX: required by collate_counts.py to build coldata
    log:
        f"{SCRATCH}/logs/08_deseq2_prep/{{cell_line}}.log",
    resources:
        runtime = 30,
    conda:
        "../../envs/environment.yaml"
    script:
        "../scripts/collate_counts.py"


rule build_read_assignment_summary:
    """
    Compile a per-sample read-assignment table across all 3 passes and
    miRNA screening. This is QC metric 2 from the proposal: proportional
    read assignment across mature tRNA / pre-tRNA / miRNA / unassigned.

    Inputs actually used:
      - Pass 1 filter stats (.filter_stats.tsv)   → pass 1 aligned / functional counts
      - Bowtie2 stats files (.bowtie2_stats.txt)  → pass 2 mapped counts
      - miRDeep2 qc_summary files                 → pass 3 miRNA read counts

    NOTE: Trim Galore trimming reports are NOT parsed here. Total input read
    counts are derived from pass1 total_aligned + unmapped (filter_stats.tsv).
    If pre-trim read counts are required, add the trimming reports as inputs
    and parse them in the run: block below.
    """
    input:
        filter_stats  = lambda wildcards: [
            f"{SCRATCH}/pass1_filters/{s}/{s}.filter_stats.tsv"
            for s in samples_for(wildcards.cell_line)
        ],
        bowtie2_stats = lambda wildcards: [
            f"{SCRATCH}/pass2_pretRNA/{s}/{s}.bowtie2_stats.txt"
            for s in samples_for(wildcards.cell_line)
        ],
        mirna_frac    = lambda wildcards: [
            f"{SCRATCH}/pass3_mirna/{s}/{s}_miRNA_fraction.tsv"
            for s in samples_for(wildcards.cell_line)
        ],
    output:
        summary = f"{SCRATCH}/qc/read_assignment/{{cell_line}}_read_assignment_summary.tsv",
    params:
        cell_line = "{cell_line}",
        scratch   = SCRATCH,
        qc_thresh = config["qc_thresholds"],
    log:
        f"{SCRATCH}/logs/08_read_assignment/{{cell_line}}.log",
    conda:
        "../../envs/environment.yaml"
    resources:
        runtime   = config["resources"]["build_read_assignment_summary"]["runtime_min"],
        sge_extra = sge_extra("build_read_assignment_summary"),
    script:
        "../scripts/build_read_assignment_summary.py"
