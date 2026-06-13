# =============================================================================
# workflow/rules/01_trim.smk
# Adapter trimming with Trim Galore (wrapping Cutadapt).
#
# KEY DESIGN: A549 has 2 sequencing lanes per sample (uncompressed .fastq).
#             THP1 has 1 lane per sample (.fq.gz).
#
#             Cutadapt ≥3.x uses dnaio's chunk-based parallel reader when
#             --cores > 1.  That reader requires seekable file descriptors.
#             Bash process substitutions (<(cat ...)) produce /dev/fd/N pipes,
#             which are NOT seekable — causing dnaio.UnknownFileFormat on all
#             multi-core runs.
#
#             Fix: detect multi-lane samples (A549) at runtime and pre-merge
#             the two lane files into a single temp file under {scratch}/tmp
#             (NOT /tmp — parallel uncompressed merges fill the small tmpfs)
#             before calling Trim Galore.  Single-lane gzipped samples (THP1)
#             are passed directly — no process substitution, no temp file.
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


rule trim_galore:
    """
    Run Trim Galore in paired-end mode.
    - Trims Illumina Universal adapter (R1) and Small RNA adapter (R2)
    - Minimum post-trim length: 15 nt (retains tRFs per proposal Section 3.2)
    - Q20 5-prime quality trimming
    - A549 (n_lanes=2): pre-merges lane files to {scratch}/tmp before trimming
      to avoid the dnaio seekability issue with process substitutions + --cores>1
      and to avoid exhausting the small /tmp tmpfs with parallel uncompressed FASTQs
    - THP1 (n_lanes=1): passes .fq.gz directly — no temp file needed
    - --basename sets output file prefix to sample_id regardless of input path
    """
    input:
        r1_files = trim_r1_input,
        r2_files = trim_r2_input,
    output:
        r1_trim   = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_val_1.fq.gz",
        r2_trim   = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_val_2.fq.gz",
        report_r1 = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_R1_trimming_report.txt",
        report_r2 = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_R2_trimming_report.txt",
    params:
        outdir     = f"{SCRATCH}/trimmed/{{sample}}",
        n_lanes    = trim_n_lanes,          # int: 1 (THP1) or 2 (A549)
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
        tmpdir    = f"{SCRATCH}/tmp",
        sge_pe    = "sharedmem",
        runtime   = 120,
        sge_extra = "-V -l h_vmem=2000M"
    conda:
        "../../envs/environment.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        mkdir -p $(dirname {log})

        echo "[$(date)] Trimming {wildcards.sample}..." > {log} 2>&1

        # ── Lane merging ────────────────────────────────────────────────────
        # Cutadapt multi-core mode uses dnaio chunk-based reading, which
        # requires seekable file descriptors.  Pre-merge multi-lane samples
        # (A549) into real temp files; pass single-lane samples (THP1) directly.
        # NOTE: temp files go to {resources.tmpdir} (on /scratch) — NOT /tmp —
        #       because parallel A549 merges of uncompressed FASTQs can exhaust
        #       the small tmpfs at /tmp (OSError: [Errno 28] No space left).
        mkdir -p {resources.tmpdir}
        if [ "{params.n_lanes}" -gt "1" ]; then
            echo "[$(date)] Merging {params.n_lanes} lanes into temp files..." >> {log}
            cat {input.r1_files} > {resources.tmpdir}/{wildcards.sample}_R1.fastq
            cat {input.r2_files} > {resources.tmpdir}/{wildcards.sample}_R2.fastq
            R1={resources.tmpdir}/{wildcards.sample}_R1.fastq
            R2={resources.tmpdir}/{wildcards.sample}_R2.fastq
        else
            R1="{input.r1_files}"
            R2="{input.r2_files}"
        fi

        # ── Trim Galore ─────────────────────────────────────────────────────
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
            "$R1" "$R2" \
            >> {log} 2>&1

        # ── Rename trimming reports to expected output names ────────────────
        # trim_galore always names reports after the INPUT filename, not
        # --basename (which only controls the trimmed read output files).
        # e.g. input A549_p4_1_R1.fastq  → A549_p4_1_R1.fastq_trimming_report.txt
        #      input sg5_c4_1_..._1.fq.gz → sg5_c4_1_..._1.fq.gz_trimming_report.txt
        # We rename to the clean {{sample}}_R1/R2_trimming_report.txt that
        # Snakemake expects.
        r1_base=$(basename "$R1")
        r2_base=$(basename "$R2")
        mv "{params.outdir}/${{r1_base}}_trimming_report.txt" "{output.report_r1}"
        mv "{params.outdir}/${{r2_base}}_trimming_report.txt" "{output.report_r2}"

        # ── Cleanup temp files (A549 only) ──────────────────────────────────
        if [ "{params.n_lanes}" -gt "1" ]; then
            rm -f {resources.tmpdir}/{wildcards.sample}_R1.fastq \
                  {resources.tmpdir}/{wildcards.sample}_R2.fastq
        fi

        echo "[$(date)] Trim Galore complete for {wildcards.sample}." >> {log}
        """
