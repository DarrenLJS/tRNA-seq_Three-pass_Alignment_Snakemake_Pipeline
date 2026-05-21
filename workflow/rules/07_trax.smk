# =============================================================================
# workflow/rules/07_trax.smk
# tRF quantification using TRAX (tRAX) from Pass 1 and Pass 2 BAMs
# (proposal Section 3.8).
#
# TRAX quantifies:
#   - 5′-tRFs (5' end fragments)
#   - 3′-tRFs (3' end fragments, including CCA-containing forms)
#   - Internal tRFs (i-tRFs)
#   - tiRNAs / tRNA halves (angiogenin-mediated cleavage at anticodon loop)
#
# tiRNA accumulation is of particular biological interest as a marker of
# stress-induced translational suppression during the antiviral response.
#
# Reads with < 3 uniquely-mapping reads across all samples are discarded
# (proposal Section 3.8).
#
# TRAX is run per cell line (all 15 samples together) for consistent tRF
# clustering and probabilistic multi-mapping assignment.
#
# NOTE: rule uses shell: (not run:) — shell is preferred here for bash loop logic
# The sample file and BAM merging are handled entirely in bash.
# =============================================================================

TRAX_CFG = config["trax"]
REF = config["references"]


rule trax_quantify:
    """
    Run TRAX tRF quantification for all samples of a cell line.

    Steps (all in bash):
      1. Merge Pass 1 (functional) + Pass 2 (pre-tRNA) BAMs per sample
         with samtools merge — TRAX needs one BAM per sample
      2. Write TRAX sample file (TSV: sample_name <TAB> bam_path)
      3. Run processsamples.py

    TODO: Verify TRAX version and exact CLI. Common variants:
        processsamples.py (older tRAX)
        trax processsamples (newer CLI wrapper)
    """
    input:
        pass1_bams = lambda wildcards: expand(
            f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.functional.bam",
            sample=samples_for(wildcards.cell_line)
        ),
        pass2_bams = lambda wildcards: expand(
            f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA.bam",
            sample=samples_for(wildcards.cell_line)
        ),
        trax_db    = REF["trax_ref_dir"],
    output:
        trf_counts   = f"{SCRATCH}/trax/{{cell_line}}/counts/tRF_counts.txt",
        tirna_counts = f"{SCRATCH}/trax/{{cell_line}}/counts/tiRNA_counts.txt",
        sample_file  = f"{SCRATCH}/trax/{{cell_line}}/trax_samples.txt",
    params:
        outdir      = f"{SCRATCH}/trax/{{cell_line}}",
        min_unique  = TRAX_CFG["min_unique_cov"],
        threads     = TRAX_CFG["threads"],
        scratch     = SCRATCH,
        # Pass sample list as a single space-separated string for bash to iterate
        samples     = lambda wildcards: " ".join(samples_for(wildcards.cell_line)),
    log:
        f"{SCRATCH}/logs/07_trax/{{cell_line}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/07_trax/{{cell_line}}.tsv",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/merged_bams
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        echo "[$(date)] Building merged BAMs for cell line: {wildcards.cell_line}"

        # Step 1 — Merge Pass 1 + Pass 2 BAMs per sample and write sample file
        > {output.sample_file}   # initialise empty sample file

        for SAMPLE in {params.samples}; do
            P1={params.scratch}/pass1_filters/$SAMPLE/$SAMPLE.functional.bam
            P2={params.scratch}/pass2_pretRNA/$SAMPLE/$SAMPLE.pretRNA.bam
            MERGED={params.outdir}/merged_bams/$SAMPLE.merged.bam

            echo "  Merging: $SAMPLE"
            samtools merge -f -@ 4 "$MERGED" "$P1" "$P2"
            samtools index "$MERGED"

            # Append to TRAX sample file: sample_name <TAB> bam_path
            printf "%s\t%s\n" "$SAMPLE" "$MERGED" >> {output.sample_file}
        done

        echo "[$(date)] Sample file written: {output.sample_file}"
        cat {output.sample_file}

        # Step 2 — Run TRAX
        echo "[$(date)] Running TRAX processsamples.py..."
        processsamples.py \
            --samplefile  {output.sample_file} \
            --database    {input.trax_db} \
            --outputdir   {params.outdir} \
            --mincoverage {params.min_unique} \
            --threads     {params.threads}

        echo "[$(date)] TRAX complete for {wildcards.cell_line}."
        """
