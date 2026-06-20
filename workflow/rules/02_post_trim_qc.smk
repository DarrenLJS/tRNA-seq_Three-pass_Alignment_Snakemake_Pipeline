# =============================================================================
# workflow/rules/02_post_trim_qc.smk
# Second-pass FastQC on trimmed FASTQ files, aggregated by MultiQC per cell line.
# Confirms successful adapter removal before alignment (proposal Section 3.2).
# =============================================================================


rule fastqc_post_trim:
    """
    Run FastQC on a single trimmed FASTQ file (R1 or R2).
    """
    input:
        fq = f"{SCRATCH}/trimmed/{{sample}}/{{sample}}_val_{{read}}.fq.gz",
    output:
        html = f"{SCRATCH}/qc/fastqc_post_trim/{{sample}}/{{sample}}_val_{{read}}_fastqc.html",
        zip  = f"{SCRATCH}/qc/fastqc_post_trim/{{sample}}/{{sample}}_val_{{read}}_fastqc.zip",
    params:
        outdir = f"{SCRATCH}/qc/fastqc_post_trim/{{sample}}",
    log:
        f"{SCRATCH}/logs/02_fastqc_post_trim/{{sample}}_{{read}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/02_fastqc_post_trim/{{sample}}_{{read}}.tsv",
    threads: 2
    resources:
        # FIX (config tidy-up): dropped sge_pe — it collapses to 1 SGE
        # slot under the EDDIE profile's --cores 1 regardless of
        # threads:, so it was already effectively requesting 1 slot.
        # FastQC also runs on a single input file per job here (not a
        # batch), so --threads 2 has no real parallel benefit to lose by
        # making that 1-slot allocation explicit instead of misleading.
        runtime   = config["resources"]["fastqc_post_trim"]["runtime_min"],
        sge_extra = sge_extra("fastqc_post_trim"),
    conda:
        "../../envs/environment.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        fastqc \
            --outdir  {params.outdir} \
            --threads {threads} \
            --quiet \
            {input.fq} \
            > {log} 2>&1
        """


def post_trim_fastqc_for(cell_line):
    """
    Collect all post-trim FastQC zip files for a given cell line
    (both R1 and R2 for every sample).
    """
    zips = []
    for s in samples_for(cell_line):
        for read in ["1", "2"]:
            zips.append(
                f"{SCRATCH}/qc/fastqc_post_trim/{s}/{s}_val_{read}_fastqc.zip"
            )
    return zips


rule multiqc_post_trim:
    """
    Aggregate per-sample post-trim FastQC reports into a single MultiQC report.
    Run per cell line so A549 and THP1 are compared separately
    (they have different read lengths: 50 bp vs 150 bp).
    """
    input:
        zips = lambda wildcards: post_trim_fastqc_for(wildcards.cell_line),
    output:
        html = f"{SCRATCH}/qc/multiqc_post_trim/{{cell_line}}_post_trim_multiqc.html",
        data = directory(
            f"{SCRATCH}/qc/multiqc_post_trim/{{cell_line}}_post_trim_multiqc_data"
        ),
    params:
        indir  = f"{SCRATCH}/qc/fastqc_post_trim",
        outdir = f"{SCRATCH}/qc/multiqc_post_trim",
        title  = "{cell_line} post-trim QC",
        fname  = "{cell_line}_post_trim_multiqc",
    log:
        f"{SCRATCH}/logs/02_multiqc_post_trim/{{cell_line}}.log",
    resources:
        runtime = 60,
    conda:
        "../../envs/environment.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})

        # Collect only the FastQC zips for this cell line into a temp filelist
        # so MultiQC doesn't pick up the other cell line's files.
        TMPLIST=$(mktemp)
        printf '%s\n' {input.zips} > "$TMPLIST"

        multiqc \
            --file-list "$TMPLIST" \
            --outdir    {params.outdir} \
            --title     "{params.title}" \
            --filename  {params.fname} \
            --force \
            > {log} 2>&1

        rm -f "$TMPLIST"
        echo "[$(date)] MultiQC post-trim done for {wildcards.cell_line}." >> {log}
        """
