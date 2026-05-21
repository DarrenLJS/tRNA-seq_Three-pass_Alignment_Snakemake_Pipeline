# =============================================================================
# workflow/rules/06_pass3_mirna.smk
# Pass 3: screen reads unmapped in both Pass 1 AND Pass 2 for miRNA content
# using miRDeep2 (proposal Section 3.3, Pass 3).
#
# Purpose: QC metric for library composition.
#   - Quantify the fraction of reads that are co-purified miRNAs
#   - Flag libraries with unexpectedly high miRNA fractions (>20%)
#     as potential co-purification anomalies (proposal Section 4, QC metric 2)
#
# IMPORTANT: miRDeep2 in quantification mode (miRDeep2_quant.pl or miRNA_quant)
# requires reads to be collapsed first (mapper.pl). The pipeline uses:
#   1. mapper.pl   — converts FASTQ to collapsed FASTA + mapping
#   2. quantifier.pl — quantifies against human miRNA hairpin/mature references
#
# For paired-end reads, miRDeep2 typically uses R1 only (the 5' read, which
# contains the miRNA seed sequence). R2 is not used here.
# =============================================================================

REF = config["references"]
MD2 = config["mirdeep2"]


rule mirdeep2_pass3:
    """
    Quantify miRNA content in Pass-2-unmapped reads using miRDeep2.

    Steps:
      1. mapper.pl  — collapse reads and map to GRCh38 with Bowtie 1
      2. quantifier.pl — quantify against hsa miRNA references

    Uses the pre-built Bowtie 1 genome index at:
      config["references"]["bowtie1_genome_index"]
    This index was built once manually — do NOT rebuild per sample.

    Only R1 reads are used (5' end carries miRNA seed sequence).
    """
    input:
        r1           = f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA_unmapped_R1.fq.gz",
        mature_hsa   = REF["mirbase_mature_hsa"],
        hairpin_hsa  = REF["mirbase_hairpin_hsa"],
    output:
        counts     = f"{SCRATCH}/pass3_mirna/{{sample}}/{{sample}}_miRNA_counts.txt",
        collapsed  = f"{SCRATCH}/pass3_mirna/{{sample}}/{{sample}}_collapsed.fa",
        qc_summary = f"{SCRATCH}/pass3_mirna/{{sample}}/{{sample}}_miRNA_fraction.tsv",
    params:
        outdir          = f"{SCRATCH}/pass3_mirna/{{sample}}",
        threads         = MD2["threads"],
        genome_bt1_idx  = REF["bowtie1_genome_index"],
        sample          = "{sample}",
    log:
        f"{SCRATCH}/logs/06_pass3_mirna/{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/06_pass3_mirna/{{sample}}.tsv",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        cd {params.outdir}

        echo "[$(date)] Pass 3 miRDeep2 for {params.sample}..." > {log}

        # Step 1: Decompress R1 (mapper.pl requires uncompressed input)
        zcat {input.r1} > {params.outdir}/{params.sample}_R1.fq

        # Step 2: mapper.pl — collapse reads and map to genome
        # -p: pre-built Bowtie 1 index prefix (do NOT rebuild)
        mapper.pl \
            {params.outdir}/{params.sample}_R1.fq \
            -e \
            -h \
            -i \
            -j \
            -m \
            -p {params.genome_bt1_idx} \
            -s {output.collapsed} \
            -t {params.outdir}/{params.sample}_mapped.arf \
            -v \
            >> {log} 2>&1

        # Step 3: quantifier.pl — quantify miRNA expression
        quantifier.pl \
            -p {input.hairpin_hsa} \
            -m {input.mature_hsa} \
            -r {output.collapsed} \
            -y {params.sample} \
            -d \
            >> {log} 2>&1

        # Move output to expected location
        if [ -f "miRNAs_expressed_{params.sample}.csv" ]; then
            mv "miRNAs_expressed_{params.sample}.csv" {output.counts}
        else
            echo "WARNING: miRDeep2 quantifier output not found." >> {log}
            echo -e "#miRNA\tread_count\tnorm_count" > {output.counts}
        fi

        # Step 4: compute miRNA fraction of total unmapped reads
        TOTAL_INPUT=$(zcat {input.r1} | awk 'NR%4==1' | wc -l)
        MIRNA_READS=$(awk -F'\t' 'NR>1 {{sum+=$2}} END {{print sum+0}}' {output.counts})
        echo -e "sample\ttotal_pass2_unmapped\tmirna_reads\tmirna_fraction" \
            > {output.qc_summary}
        awk -v s="{params.sample}" \
            -v t="$TOTAL_INPUT" \
            -v m="$MIRNA_READS" \
            'BEGIN {{
                frac = (t > 0) ? m/t : 0;
                printf "%s\t%d\t%d\t%.4f\n", s, t, m, frac
            }}' >> {output.qc_summary}

        # Clean up temporary uncompressed file
        rm -f {params.outdir}/{params.sample}_R1.fq

        echo "[$(date)] Pass 3 miRDeep2 complete for {params.sample}." >> {log}
        """
