# =============================================================================
# workflow/rules/08_data_prep.smk
# Collate mim-tRNAseq count tables into per-cell-line count matrices and
# write the corresponding coldata (sample metadata) TSVs for DESeq2/edgeR.
# Also assembles wobble-position mismatch matrices for the Binomial GLM
# (future rule 09) and a per-cell-line read assignment summary (QC).
#
# Outputs per cell line:
#   isodecoder_counts_matrix.tsv   — rows=isodecoders, cols=samples
#   isoacceptor_counts_matrix.tsv  — rows=isoacceptors, cols=samples
#   coldata.tsv                    — sample metadata for DESeq2
#   pos34_coverage_matrix.tsv      — wobble-pos coverage (GLM trials)
#   pos34_mismatch_matrix.tsv      — wobble-pos mismatch count (GLM successes)
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
        # Per-position misincorporation data from mim-tRNAseq (rule 03).
        # Used to compute wobble-position 34 coverage QC columns:
        #   n_isodecoders_pos34_covered, median_pos34_coverage, flag_low_mismatch_cov.
        mismatch_dir  = lambda wildcards: (
            f"{SCRATCH}/pass1_mimtrnaseq/{wildcards.cell_line}/mismatch"
        ),
    output:
        summary = f"{SCRATCH}/qc/read_assignment/{{cell_line}}_read_assignment_summary.tsv",
    params:
        cell_line     = "{cell_line}",
        scratch       = SCRATCH,
        qc_thresh     = config["qc_thresholds"],
        min_pos34_cov = config["qc_thresholds"]["min_pos34_coverage"],
    log:
        f"{SCRATCH}/logs/08_read_assignment/{{cell_line}}.log",
    conda:
        "../../envs/environment.yaml"
    resources:
        runtime   = config["resources"]["build_read_assignment_summary"]["runtime_min"],
        sge_extra = sge_extra("build_read_assignment_summary"),
    script:
        "../scripts/build_read_assignment_summary.py"


rule build_mismatch_matrices:
    """
    Assemble per-cell-line wobble-position (34) mismatch matrices for
    the Binomial GLM wobble-modification inference (future rule 09).

    Directly analogous to build_deseq2_inputs for count matrices:
        mismatch/ → pos34_coverage_matrix.tsv  (trials: read depth at pos 34)
                  → pos34_mismatch_matrix.tsv   (successes: mismatch count at pos 34)

    Both matrices: rows = isodecoder clusters, cols = samples (manifest order).
    If mim-tRNAseq stores proportions rather than raw counts, collate_mismatch.py
    recovers counts as round(proportion * coverage).

    Calls workflow/scripts/collate_mismatch.py.
    """
    input:
        mismatch_dir = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/mismatch",
    output:
        cov_matrix = f"{SCRATCH}/deseq2_input/{{cell_line}}/pos34_coverage_matrix.tsv",
        mm_matrix  = f"{SCRATCH}/deseq2_input/{{cell_line}}/pos34_mismatch_matrix.tsv",
    params:
        outdir    = f"{SCRATCH}/deseq2_input/{{cell_line}}",
        cell_line = "{cell_line}",
        manifest  = config["manifest"],
    log:
        f"{SCRATCH}/logs/08_mismatch_matrices/{{cell_line}}.log",
    resources:
        runtime = 30,
    conda:
        "../../envs/environment.yaml"
    script:
        "../scripts/collate_mismatch.py"
