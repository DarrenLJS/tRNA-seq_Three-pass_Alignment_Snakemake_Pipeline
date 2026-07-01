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

MD2 = config["mirdeep2"]
# NOTE: REF is defined in 00_reference_prep.smk (same included namespace).
#       Accessing config["references"] directly to avoid redefinition.


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

    If zero reads survive pass2 (all reads mapped to pre-tRNA), the rule
    exits cleanly with empty outputs rather than crashing in mapper.pl.
    """
    input:
        r1           = f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA_unmapped_R1.fq.gz",
        mature_hsa   = config["references"]["mirbase_mature_hsa"],
        hairpin_hsa  = config["references"]["mirbase_hairpin_hsa"],
    output:
        counts     = f"{SCRATCH}/pass3_mirna/{{sample}}/{{sample}}_miRNA_counts.txt",
        collapsed  = f"{SCRATCH}/pass3_mirna/{{sample}}/{{sample}}_collapsed.fa",
        qc_summary = f"{SCRATCH}/pass3_mirna/{{sample}}/{{sample}}_miRNA_fraction.tsv",
    params:
        outdir          = f"{SCRATCH}/pass3_mirna/{{sample}}",
        threads         = MD2["threads"],
        genome_bt1_idx  = config["references"]["bowtie1_genome_index"],
        sample          = "{sample}",
    log:
        f"{SCRATCH}/logs/06_pass3_mirna/{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/06_pass3_mirna/{{sample}}.tsv",
    threads: lambda wildcards: MD2["threads"]
    resources:
        # FIX (config tidy-up): dropped sge_pe — it collapses to 1 SGE
        # slot under the EDDIE profile's --cores 1 regardless of
        # threads:, which was already the correct allocation here since
        # mapper.pl and quantifier.pl are effectively single-threaded
        # (mirdeep2.threads is unused by either tool). Memory was the
        # real fix needed (was 2000M — Bowtie1 mapping against hg38
        # peaked at 2.54G vmem, exceeding the 2G limit and triggering an
        # SGE OOM kill, exit 137, failed 46 — raised to 6000M).
        # Vmem/runtime now live in config["resources"]["mirdeep2_pass3"].
        runtime   = config["resources"]["mirdeep2_pass3"]["runtime_min"],
        sge_extra = sge_extra("mirdeep2_pass3"),
    conda:
        "../../envs/environment.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir} $(dirname {log})
        cd {params.outdir}

        echo "[$(date)] Pass 3 miRDeep2 for {params.sample}..." > {log}

        # FIX: guard against empty pass2 unmapped output.
        # mapper.pl crashes with uninitialised-value Perl warnings when fed an
        # empty FASTQ (all reads were absorbed by pass2 pre-tRNA alignment).
        # Write valid empty outputs and exit 0 so downstream rules still run.
        READ_COUNT=$(zcat {input.r1} | awk 'NR%4==1' | wc -l)
        if [ "$READ_COUNT" -eq 0 ]; then
            echo "[$(date)] Zero reads in pass2 unmapped input — skipping miRDeep2 for {params.sample}." >> {log}
            touch {output.collapsed}
            echo -e "#miRNA\tread_count\tnorm_count" > {output.counts}
            echo -e "sample\ttotal_pass2_unmapped\tmirna_reads\tmirna_fraction" > {output.qc_summary}
            echo -e "{params.sample}\t0\t0\t0.0000" >> {output.qc_summary}
            exit 0
        fi

        # Step 1: Decompress R1 (mapper.pl requires uncompressed input)
        zcat {input.r1} > {params.outdir}/{params.sample}_R1.fq

        # Step 2: mapper.pl — collapse reads and map to genome
        # -p: pre-built Bowtie 1 index prefix (do NOT rebuild)
        # FIX: mapper.pl refuses to overwrite an existing .arf file and exits
        # non-zero, causing the whole job to fail under set -euo pipefail.
        # Remove any stale .arf from a previous partial run before proceeding.
        rm -f {params.outdir}/{params.sample}_mapped.arf
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
        # FIX: quantifier.pl (invoked with -y {params.sample}) names its real
        # output "miRNAs_expressed_all_samples_{sample}.csv", NOT
        # "miRNAs_expressed_{sample}.csv". The old check here always missed,
        # silently falling into the else branch and writing a header-only
        # counts file for every sample -- masking substantial real miRNA
        # signal (e.g. A549_c2_1 has 1,027 expressed miRNAs / ~13.9M reads
        # in the real output, all zeroed out downstream by this bug).
        if [ -f "miRNAs_expressed_all_samples_{params.sample}.csv" ]; then
            mv "miRNAs_expressed_all_samples_{params.sample}.csv" {output.counts}
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

        # Clean up temporary uncompressed file and mapper.pl working directory.
        # mapper.pl creates a dir_mapper_seq_*/ temp dir alongside the input fq
        # and does not remove it on exit. On a quota-limited filesystem these
        # accumulate across samples and cause subsequent runs to fail with
        # "Disk quota exceeded".
        rm -f {params.outdir}/{params.sample}_R1.fq
        rm -rf {params.outdir}/dir_mapper_seq_*

        echo "[$(date)] Pass 3 miRDeep2 complete for {params.sample}." >> {log}
        """
