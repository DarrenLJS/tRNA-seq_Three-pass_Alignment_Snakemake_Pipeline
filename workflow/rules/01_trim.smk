# =============================================================================
# workflow/rules/01_trim.smk
# Adapter trimming with Trim Galore (wrapping Cutadapt).
#
# KEY DESIGN: A549 has 2 sequencing lanes per sample (uncompressed .fastq).
#             THP1 has 1 lane per sample (.fq.gz).
#
#             Rather than pre-merging A549 lanes (slow, ~5 hr for one sample),
#             we stream them on-the-fly using bash process substitution:
#               <(cat lane1_R1.fastq lane2_R1.fastq)
#             Trim Galore sees a single continuous stream; no disk I/O overhead.
#             The --basename flag names output files by sample_id regardless
#             of the input path (critical when using /dev/fd/N from process sub).
#
# Outputs (for every sample):
#   {scratch}/trimmed/{sample}/{sample}_val_1.fq.gz   (trimmed R1)
#   {scratch}/trimmed/{sample}/{sample}_val_2.fq.gz   (trimmed R2)
#   {scratch}/trimmed/{sample}/{sample}_R1_trimming_report.txt
#   {scratch}/trimmed/{sample}/{sample}_R2_trimming_report.txt
# =============================================================================

TG = config["trim_galore"]


def trim_r1_input(wildcards):
    """Return list of raw R1 FASTQ paths — used as Snakemake input for DAG."""
    return r_files(wildcards.sample, "R1")

def trim_r2_input(wildcards):
    """Return list of raw R2 FASTQ paths."""
    return r_files(wildcards.sample, "R2")

def trim_r1_cmd(wildcards):
    """Process-substitution expression for R1 (handles single or multi-lane)."""
    return cat_cmd(wildcards.sample, "R1")

def trim_r2_cmd(wildcards):
    """Process-substitution expression for R2."""
    return cat_cmd(wildcards.sample, "R2")


rule trim_galore:
    """
    Run Trim Galore in paired-end mode.
    - Trims Illumina Universal adapter (R1) and Small RNA adapter (R2)
    - Minimum post-trim length: 15 nt (retains tRFs per proposal Section 3.2)
    - Q20 5-prime quality trimming
    - A549: streams 2 lane files per read via <(cat ...)  — no disk merge
    - THP1: streams single .fq.gz via <(zcat ...)
    - --basename sets output file prefix to sample_id
    """
    input:
        r1_files = trim_r1_input,
        r2_files = trim_r2_input,
    output:
        r1_trim  = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_val_1.fq.gz",
        r2_trim  = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_val_2.fq.gz",
        report_r1 = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_R1_trimming_report.txt",
        report_r2 = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_R2_trimming_report.txt",
    params:
        outdir     = f"{SCRATCH}/trimmed/{{sample}}",
        r1_cmd     = trim_r1_cmd,
        r2_cmd     = trim_r2_cmd,
        adapter_r1 = TG["adapter_r1"],
        adapter_r2 = TG["adapter_r2"],
        min_len    = TG["min_length"],
        quality    = TG["quality"],
        cores      = TG["cores"],
        basename   = "{sample}",
    log:
        f"{SCRATCH}/logs/01_trim/{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/01_trim/{{sample}}.tsv",
    threads: lambda wildcards: TG["cores"]
    resources:
        mem_mb = config["resources"]["trim_mem_mb"],
    shell:
        # NOTE: process substitution requires bash (Snakemake uses bash by default)
        # The --basename flag (Trim Galore ≥0.6.0) is required here so that
        # output files are named {sample}_val_1/2.fq.gz regardless of the
        # /dev/fd/N path that process substitution creates.
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        echo "[$(date)] Trimming {wildcards.sample}..." > {log} 2>&1

        trim_galore \
            --paired \
            --adapter   {params.adapter_r1} \
            --adapter2  {params.adapter_r2} \
            --length    {params.min_len} \
            --quality   {params.quality} \
            --cores     {params.cores} \
            --gzip \
            --basename  {params.basename} \
            --output_dir {params.outdir} \
            {params.r1_cmd} \
            {params.r2_cmd} \
            >> {log} 2>&1

        echo "[$(date)] Trim Galore complete for {wildcards.sample}." >> {log}
        """
