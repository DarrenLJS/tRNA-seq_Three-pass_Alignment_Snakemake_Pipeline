# =============================================================================
# workflow/rules/05_pass2_pretRNA.smk
# Pass 2: align reads unmapped in Pass 1 (from mim-tRNAseq) to the pre-tRNA
# reference using Bowtie2 in end-to-end sensitive mode (proposal Section 3.3).
#
# The pre-tRNA reference contains:
#   - Intron-retaining forms (for reads spanning intron–exon junctions)
#   - Spliced (exon-only) forms (for mature exon reads that escaped Pass 1)
#
# Unmapped reads from this pass are forwarded to Pass 3 (miRDeep2).
#
# Key outputs per sample:
#   {sample}.pretRNA.bam          — aligned reads (pre-tRNA mapped)
#   {sample}.pretRNA_unmapped_R1.fq.gz  → input to Pass 3
#   {sample}.pretRNA_unmapped_R2.fq.gz  → input to Pass 3
# =============================================================================

BT2 = config["bowtie2"]
REF = config["references"]


rule bowtie2_pretRNA:
    """
    Align Pass-1-unmapped reads to the pre-tRNA reference with Bowtie2.
    Uses --end-to-end --sensitive mode and allows up to 20 alignments
    (-k 20) to handle the redundancy inherent in pre-tRNA sequences.
    Outputs: sorted+indexed BAM of pre-tRNA alignments, plus unmapped FASTQ pair.
    """
    input:
        r1          = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.unmapped_R1.fq.gz",
        r2          = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.unmapped_R2.fq.gz",
        index_done  = REF["pretRNA_index"] + ".1.bt2",
    output:
        bam          = f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA.bam",
        bai          = f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA.bam.bai",
        unmapped_r1  = f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA_unmapped_R1.fq.gz",
        unmapped_r2  = f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA_unmapped_R2.fq.gz",
        align_stats  = f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.bowtie2_stats.txt",
    params:
        outdir        = f"{SCRATCH}/pass2_pretRNA/{{sample}}",
        index_prefix  = REF["pretRNA_index"],
        preset        = BT2["preset"],
        max_align     = BT2["max_align"],
        threads       = BT2["threads"],
    log:
        f"{SCRATCH}/logs/05_pass2_pretRNA/{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/05_pass2_pretRNA/{{sample}}.tsv",
    threads: lambda wildcards: BT2["threads"]
    resources:
        mem_mb = config["resources"]["bowtie2_mem_mb"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})

        echo "[$(date)] Bowtie2 Pass 2 for {wildcards.sample}..." > {log}

        # Run Bowtie2; pipe through samtools to sort and extract unmapped pairs
        bowtie2 \
            {params.preset} \
            -k {params.max_align} \
            -p {params.threads} \
            --no-mixed \
            --no-discordant \
            -x {params.index_prefix} \
            -1 {input.r1} \
            -2 {input.r2} \
            --un-conc-gz {params.outdir}/{wildcards.sample}.pretRNA_unmapped_R%.fq.gz \
            2> {output.align_stats} \
        | samtools sort -@ 4 -o {output.bam} -

        # Rename --un-conc-gz outputs to match expected naming convention
        # Bowtie2 --un-conc-gz with % wildcard outputs _R1.fq.gz and _R2.fq.gz
        # (already matching our expected output names above)

        samtools index {output.bam}

        echo "[$(date)] Bowtie2 Pass 2 done for {wildcards.sample}." >> {log}
        cat {output.align_stats} >> {log}
        """
