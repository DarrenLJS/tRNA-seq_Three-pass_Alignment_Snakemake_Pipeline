# =============================================================================
# workflow/rules/07_trax.smk
# tRF quantification using TRAX (tRAX) from Pass 1 and Pass 2 BAMs
# (proposal Section 3.8).
#
# TRAX quantifies:
#   - 5'-tRFs (5' end fragments)
#   - 3'-tRFs (3' end fragments, including CCA-containing forms)
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
# NOTE: rule uses shell: (not run:) -- shell is preferred here for bash loop
# logic. The sample file and BAM merging are handled entirely in bash.
#
# FIX HISTORY (2026-06-13):
#   1. h_vmem raised 4000M -> 16000M (OOM kills on busy nodes)
#   2. samplefile format: 2-col -> 3-col (name, condition, bam_path) required
#      by trnasequtils.samplefile() which reads fields[2] for the BAM path
#   3. processsamples.py arguments corrected:
#        --database  -> --databasename
#        --threads   -> --cores
#        --outputdir -> removed (tRAX has no --outputdir; outputs go to CWD)
#        (new)       -> --experimentname (required)
#        (new)       -> cd to outdir before running so outputs land there
# FIX HISTORY (2026-06-16):
#   8. Wall-clock timeout: both cell lines hit h_rt=24h mid-bowtie2 (~18.5 h
#      elapsed, mapping still in progress for THP1). runtime raised from
#      1440 -> 4320 min (72 h) and pulled from
#      config["resources"]["trax_runtime_min"] so it can be adjusted without
#      editing this file.
# FIX HISTORY (2026-06-15):
#   7. Disk quota failures on rerun:
#      Merged BAMs + FASTQs (~100 GB) were never deleted after a successful
#      tRAX run. On the next Snakemake rerun the merge loop re-created them on
#      top of the existing files, exhausting disk quota and causing samtools to
#      report the misleading error "failed writing: No such file or directory".
#      Fix A: fast-path check moved to BEFORE the merge loop so a rerun skips
#      all BAM/FASTQ I/O entirely when outputs already exist.
#      Fix B: merged BAMs, index files, and FASTQs are deleted after a
#      successful tRAX run since they are pure intermediates.
# FIX HISTORY (2026-06-14):
#   4. SGE RSS kill (exit 137, failed 47: execd enforced h_rss limit):
#      Snakemake head-node runs with --cores 1, which caps job.threads to 1
#      for SGE slot submission regardless of the rule's threads: value.
#      Result: SGE submitted -pe sharedmem 1 -l h_vmem=32000M (32G for 1 slot)
#      but tRAX spawned 8 parallel bowtie2 processes whose combined RSS (~17G)
#      exhausted the single-slot RSS budget and was killed.
#      Fix: remove sge_pe (which causes the profile to emit -pe sge_pe threads
#      and therefore inherits the capped slot count) and embed the full PE
#      specification directly in sge_extra as -pe sharedmem 8.
#      8 slots x 8000M h_vmem = 64G total RSS budget -- sufficient headroom.
#   5. bowtie2 "reads file does not look like a FASTQ file" / SIGABRT:
#      tRAX calls bowtie2 with -U <file>, which expects FASTQ input, not an
#      aligned BAM.  We were passing merged aligned BAMs directly to tRAX,
#      which caused bowtie2 to SIGABRT on every sample.
#      Fix: convert each merged BAM to FASTQ with samtools fastq before tRAX.
#   6. Wrong FASTQ extraction -- all pass1 reads silently discarded:
#      The merged BAM contains TWO DIFFERENT read types:
#        - Pass 1 (functional): single-end reads, NO pairing flags set.
#          mim-tRNAseq takes R1 only and produces .unpaired_uniq.bam;
#          the reads have no 0x1/0x40/0x80 flags.
#        - Pass 2 (pre-tRNA): paired-end reads, proper 0x40/0x80 flags set.
#          bowtie2 was run with -1 R1 -2 R2 --no-mixed --no-discordant.
#      The previous extraction used samtools fastq -1 R1 -2 R2 -0 /dev/null,
#      which routes reads with no pairing flag to /dev/null -- silently
#      discarding the entire pass1 functional set (the majority of reads).
#      This caused bowtie2 to see empty or near-empty FASTQ files and report
#      "Unable to read file magic number / 0 reads / 0% alignment rate".
#      Fix: use plain "samtools fastq MERGED > READS" which outputs all reads
#      unconditionally, regardless of pairing flags. Since tRAX maps each read
#      independently via bowtie2 -U, no paired-end handling is needed here.
#      The name-sort step is also removed as it served only the old -1/-2 mode.
# =============================================================================

TRAX_CFG = config["trax"]
# NOTE: REF is defined in 00_reference_prep.smk (same included namespace).
#       Accessing config["references"] directly to avoid redefinition.


rule trax_quantify:
    """
    Run TRAX tRF quantification for all samples of a cell line.

    Steps (all in bash):
      1. Merge Pass 1 (functional) + Pass 2 (pre-tRNA) BAMs per sample.
         Pass 1 BAMs are single-end (mim-tRNAseq, R1 only, no pairing flags).
         Pass 2 BAMs are paired-end (bowtie2 -1/-2, flags 0x40/0x80 set).
      2. Convert the merged BAM to FASTQ unconditionally with
         "samtools fastq MERGED > READS" -- all reads output regardless of
         pairing flags. tRAX maps each read independently via bowtie2 -U,
         so no paired-end handling is needed.
      3. Write TRAX sample file (TSV: name TAB condition TAB fastq_path)
      4. Run processsamples.py

    tRAX scripts are called via:  python {params.trax_dir}/processsamples.py
    because tRAX is a GitHub repo, not a conda/pip package.  Its dependencies
    (samtools, pysam, bowtie2) are provided by trax_env.yaml.
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
        trax_db      = config["references"]["trax_ref_dir"],
        ensembl_gtf  = config["references"]["ensembl_gtf"],
    output:
        # FIX: corrected filenames to match actual tRAX output (verified from
        # ls trax/A549/A549/ after a completed run). The originally declared
        # names (-readcounts.txt, -trnacounts.txt, -anticodoncounts.txt) do not
        # exist; tRAX produces these instead:
        readcounts      = f"{SCRATCH}/trax/{{cell_line}}/{{cell_line}}/{{cell_line}}-normalizedreadcounts.txt",
        trnacounts      = f"{SCRATCH}/trax/{{cell_line}}/{{cell_line}}/{{cell_line}}-trnamapinfo.txt",
        anticodoncounts = f"{SCRATCH}/trax/{{cell_line}}/{{cell_line}}/{{cell_line}}-aminocounts.txt",
        sample_file  = f"{SCRATCH}/trax/{{cell_line}}/trax_samples.txt",
    params:
        outdir      = f"{SCRATCH}/trax/{{cell_line}}",
        min_unique  = TRAX_CFG["min_unique_cov"],
        threads     = TRAX_CFG["threads"],
        scratch     = SCRATCH,
        trax_dir    = config["trax"]["script_dir"],
        # Pass sample list as a single space-separated string for bash to iterate
        samples     = lambda wildcards: " ".join(samples_for(wildcards.cell_line)),
    log:
        f"{SCRATCH}/logs/07_trax/{{cell_line}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/07_trax/{{cell_line}}.tsv",
    threads: lambda wildcards: TRAX_CFG["threads"]
    resources:
        runtime   = config["resources"]["trax_runtime_min"],  # FIX 8: 1440->4320 min (72h); pulled from config
        # FIX 4: sge_pe removed -- the eddie profile emits -pe sge_pe threads,
        # but Snakemake caps threads to 1 (head node uses --cores 1), so only
        # 1 slot was granted while tRAX ran 8 bowtie2 processes -> RSS kill.
        # Hard-wiring -pe sharedmem 8 in sge_extra bypasses thread-scaling.
        # 8 slots x 8000M = 64G total RSS budget.
        sge_extra = "-V -pe sharedmem 8 -l h_vmem=8000M"
    conda:
        "../../envs/trax_env.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/merged_bams
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        echo "[$(date)] Building merged BAMs for cell line: {wildcards.cell_line}"

        # FIX 7A: fast-path moved to BEFORE the merge loop.
        # If the primary tRAX output already exists, skip all BAM/FASTQ I/O
        # and processsamples.py. Just regenerate the sample file (which Snakemake
        # deleted as a "corrupted" output) and exit 0.
        if [ -f "{output.readcounts}" ]; then
            echo "[$(date)] tRAX outputs already present — fast-path: regenerating sample file only."
            > {output.sample_file}
            for SAMPLE in {params.samples}; do
                READS={params.outdir}/merged_bams/$SAMPLE.fastq
                COND_CHAR=$(echo "$SAMPLE" | cut -d_ -f2 | cut -c1)
                COND="control"
                [ "$COND_CHAR" = "p" ] && COND="polyi"
                printf "%s\t%s\t%s\n" "$SAMPLE" "$COND" "$READS" >> {output.sample_file}
            done
            echo "[$(date)] TRAX fast-path complete for {wildcards.cell_line}."
            exit 0
        fi

        # Step 1: merge BAMs, convert to FASTQ, write sample file.
        # FIX 5: tRAX needs FASTQ not aligned BAM (bowtie2 uses -U flag).
        # FIX 6: the merged BAM is a mix of single-end reads (pass1, mim-tRNAseq
        # R1-only, no pairing flags) and paired-end reads (pass2, bowtie2 -1/-2,
        # flags 0x40/0x80 set). The previous -1/-2/-0 /dev/null extraction
        # silently dropped all pass1 reads to /dev/null because they lack pairing
        # flags. Fix: use plain samtools fastq without -1/-2 flags so all reads
        # are output unconditionally. No name-sort needed -- tRAX maps each read
        # independently via bowtie2 -U so pair order is irrelevant.
        > {output.sample_file}

        for SAMPLE in {params.samples}; do
            P1={params.scratch}/pass1_filters/$SAMPLE/$SAMPLE.functional.bam
            P2={params.scratch}/pass2_pretRNA/$SAMPLE/$SAMPLE.pretRNA.bam
            MERGED={params.outdir}/merged_bams/$SAMPLE.merged.bam
            READS={params.outdir}/merged_bams/$SAMPLE.fastq

            echo "  Merging: $SAMPLE"
            samtools merge -f -@ 4 "$MERGED" "$P1" "$P2"
            samtools index "$MERGED"

            echo "  Converting to FASTQ: $SAMPLE"
            samtools fastq -@ 4 "$MERGED" > "$READS"

            # Derive condition from sample name:
            #   _c_ prefix -> control   (e.g. A549_c2_1, THP1_c4_2)
            #   _p_ prefix -> polyi     (e.g. A549_p2_1, THP1_p8_3)
            COND_CHAR=$(echo "$SAMPLE" | cut -d_ -f2 | cut -c1)
            COND="control"
            [ "$COND_CHAR" = "p" ] && COND="polyi"

            printf "%s\t%s\t%s\n" "$SAMPLE" "$COND" "$READS" >> {output.sample_file}
        done

        echo "[$(date)] Sample file written: {output.sample_file}"
        cat {output.sample_file}

        # Step 2: run tRAX.
        # Remove per-sample BAMs left by any previous partial tRAX run.
        # tRAX checks existing BAMs in CWD against the sample-file FASTQs
        # and aborts if the counts differ. merged_bams/ subdir is unaffected.
        rm -f {params.outdir}/*.bam

        echo "[$(date)] Running TRAX processsamples.py..."
        cd {params.outdir}
        python {params.trax_dir}/processsamples.py \
            --experimentname {wildcards.cell_line} \
            --databasename   {input.trax_db}/hsapi38 \
            --samplefile     {output.sample_file} \
            --mincoverage    {params.min_unique} \
            --cores          {params.threads} \
            --ensembl        {input.ensembl_gtf}

        echo "[$(date)] TRAX complete for {wildcards.cell_line}."

        # FIX 7B: delete large intermediates after a successful tRAX run.
        # Merged BAMs + FASTQs are ~100 GB+ combined and are not needed once
        # tRAX has finished. Leaving them caused disk quota failures on reruns.
        echo "[$(date)] Cleaning up intermediate merged BAMs and FASTQs..."
        rm -f {params.outdir}/merged_bams/*.merged.bam
        rm -f {params.outdir}/merged_bams/*.merged.bam.bai
        rm -f {params.outdir}/merged_bams/*.fastq
        echo "[$(date)] Cleanup complete."
        """
