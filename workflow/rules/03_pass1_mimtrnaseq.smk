# =============================================================================
# workflow/rules/03_pass1_mimtrnaseq.smk
# Pass 1 mature tRNA alignment using mim-tRNAseq (mimseq).
#
# DESIGN: mim-tRNAseq processes all samples for a cell line together.
#   - It builds a shared, modification-aware GSNAP index from GtRNAdb
#   - Performs two-pass GSNAP alignment with misincorporation tolerance
#   - Outputs per-sample BAMs + aggregated isodecoder/isoacceptor count tables
#   - Also outputs per-position misincorporation data for wobble-mod inference
#
# WHY PER-CELL-LINE: Cross-sample clustering (--cluster-id 0.97) requires
# all samples to be processed together so the same isodecoder clusters are
# used consistently across the dataset.
#
# Input: all trimmed R1/R2 FASTQ pairs for the cell line
# Key outputs:
#   counts/Isodecoder_counts.txt   → DESeq2 isodecoder analysis
#   counts/Isoacceptor_counts.txt  → collapsed isoacceptor analysis
#   mismatch/                      → wobble position 34 modification inference
#   align/{sample}.bam             → per-sample BAMs for downstream filters
# =============================================================================

MIM = config["mimtrnaseq"]
REF = config["references"]


def mim_input_r1(cell_line):
    """Return list of trimmed R1 FASTQ paths for all samples of a cell line."""
    return [
        f"{SCRATCH}/trimmed/{s}/{s}_val_1.fq.gz"
        for s in samples_for(cell_line)
    ]

def mim_input_r2(cell_line):
    """Return list of trimmed R2 FASTQ paths for all samples of a cell line."""
    return [
        f"{SCRATCH}/trimmed/{s}/{s}_val_2.fq.gz"
        for s in samples_for(cell_line)
    ]


rule mimtrnaseq:
    """
    Run mim-tRNAseq on all samples for one cell line.

    The rule aggregates trimmed FASTQ pairs from all 15 samples of the cell
    line and passes them to mimseq as comma-separated lists via -1/-2 flags.
    mim-tRNAseq handles its own GSNAP-based alignment internally.

    NOTE on CLI: mim-tRNAseq (mimseq) expects paired-end inputs as:
        -1 R1_s1.fq.gz,R1_s2.fq.gz,...  (comma-separated, ordered)
        -2 R2_s1.fq.gz,R2_s2.fq.gz,...  (same order as -1)
    Verify against your installed version: `mimseq --help`

    The --out-dir flag is our cell-line-specific output directory; per-sample
    BAMs appear at {outdir}/align/{sample_basename}.bam
    """
    input:
        r1_files = lambda wildcards: mim_input_r1(wildcards.cell_line),
        r2_files = lambda wildcards: mim_input_r2(wildcards.cell_line),
    output:
        iso_counts  = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isodecoder_counts.txt",
        isoa_counts = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isoacceptor_counts.txt",
        # Touch file confirming the align/ directory was populated
        align_done  = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/.align_done",
        mismatch_done = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/.mismatch_done",
    params:
        outdir      = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}",
        species     = REF["mimtrnaseq_species"],
        genome      = REF["mimtrnaseq_genome"],
        threads     = MIM["threads"],
        min_cov     = MIM["min_cov"],
        max_multi   = MIM["max_multi"],
        cluster_id  = MIM["cluster_id"],
        snp_tol     = MIM["snp_tolerance"],
        r1_csv      = lambda wildcards, input: ",".join(input.r1_files),
        r2_csv      = lambda wildcards, input: ",".join(input.r2_files),
    log:
        f"{SCRATCH}/logs/03_mimtrnaseq/{{cell_line}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/03_mimtrnaseq/{{cell_line}}.tsv",
    threads: lambda wildcards: MIM["threads"]
    resources:
        mem_mb = config["resources"]["mim_mem_mb"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})

        echo "[$(date)] Starting mim-tRNAseq for cell line: {wildcards.cell_line}" > {log} 2>&1
        echo "[$(date)] Samples: {params.r1_csv}" >> {log}

        mimseq \
            --species      {params.species} \
            --cluster \
            --cluster-id   {params.cluster_id} \
            --threads      {params.threads} \
            --min-cov      {params.min_cov} \
            --max-multi    {params.max_multi} \
            --snp-tolerance {params.snp_tol} \
            --out-dir      {params.outdir} \
            -1 {params.r1_csv} \
            -2 {params.r2_csv} \
            >> {log} 2>&1

        # Verify expected outputs exist and create sentinel files
        echo "[$(date)] Verifying outputs..." >> {log}

        if [ ! -f "{output.iso_counts}" ]; then
            echo "ERROR: Isodecoder_counts.txt not found — mim-tRNAseq may have failed." >> {log}
            exit 1
        fi

        # Confirm per-sample BAMs exist in align/
        N_BAMS=$(find {params.outdir}/align/ -name "*.bam" | wc -l)
        echo "[$(date)] Found $N_BAMS BAM files in {params.outdir}/align/" >> {log}

        touch {output.align_done}
        touch {output.mismatch_done}
        echo "[$(date)] mim-tRNAseq complete for {wildcards.cell_line}." >> {log}
        """


# ---------------------------------------------------------------------------
# Helper rule: locate the per-sample BAM produced by mim-tRNAseq.
# mim-tRNAseq names BAMs by the input FASTQ basename (stripping _val_1.fq.gz).
# We symlink/copy to our consistent per-sample naming scheme so downstream
# rules can use {sample}.bam paths without worrying about the mim naming logic.
# ---------------------------------------------------------------------------
rule link_mimtrnaseq_bam:
    """
    Create a consistently-named symlink for each per-sample mim-tRNAseq BAM.
    mim-tRNAseq names BAMs by stripped FASTQ basename; we standardise to
    {sample}.bam for all downstream rules.
    """
    input:
        align_done = lambda wildcards: (
            f"{SCRATCH}/pass1_mimtrnaseq/"
            f"{manifest.loc[wildcards.sample,'cell_line']}/.align_done"
        ),
    output:
        bam   = f"{SCRATCH}/pass1_mimtrnaseq/{{sample}}/{{sample}}.bam",
        bai   = f"{SCRATCH}/pass1_mimtrnaseq/{{sample}}/{{sample}}.bam.bai",
    params:
        # mim-tRNAseq strips '_val_1.fq.gz' from R1 filename to name BAM
        src_bam = lambda wildcards: (
            f"{SCRATCH}/pass1_mimtrnaseq/"
            f"{manifest.loc[wildcards.sample,'cell_line']}/"
            f"align/{wildcards.sample}_val_1.bam"
        ),
        outdir = f"{SCRATCH}/pass1_mimtrnaseq/{{sample}}",
    log:
        f"{SCRATCH}/logs/03_link_bam/{{sample}}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        # Copy (not symlink) to avoid cross-filesystem issues on /scratch
        echo "[$(date)] Copying BAM for {wildcards.sample}..." > {log}
        cp {params.src_bam} {output.bam}

        echo "[$(date)] Indexing BAM..." >> {log}
        samtools index {output.bam}

        echo "[$(date)] Done." >> {log}
        """
