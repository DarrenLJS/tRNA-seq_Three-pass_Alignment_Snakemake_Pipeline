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
# FIX HISTORY (2026-06-18):
#   9. Mapping phase ran in two sequential waves: TRAX's --cores controls how
#      many samples map concurrently, but each internal bowtie2 call is
#      single-threaded (-p 1) -- TRAX has no per-sample multithreading lever.
#      With --cores 8 and 15 samples/cell line, mapping ran as 8-then-7
#      sequential waves, contributing the bulk of a 22.5h mapping+counting
#      phase (A549: 09:40 start -> next-day 08:19 complete). Raised
#      trax.threads 8->16 (config) and -pe sharedmem 8->16 (this file) so all
#      15 samples map in a single wave. h_vmem kept at 8000M/slot (128G total
#      budget); prior peak was ~44G against a 64G budget at 8 slots, so this
#      preserves the same headroom ratio at 16. Eddie nodes are 32-128 core
#      (confirmed via qhost), so 16 sharedmem slots carries no queue risk.
# FIX HISTORY (2026-06-16):
#   8. Wall-clock timeout: both cell lines hit h_rt=24h mid-bowtie2 (~18.5 h
#      elapsed, mapping still in progress for THP1). runtime raised from
#      1440 -> 4320 min (72 h) and pulled from
#      config["resources"]["trax_quantify"]["runtime_min"] so it can be
#      adjusted without editing this file.
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
# FIX HISTORY (2026-06-29):
#  10. THP1 TRAX job timed out three times at 72h. The p8 samples have
#      170-300M reads each; at bowtie2 -p 1 (hardcoded by TRAX) with
#      --very-sensitive -k 100, a single p8 sample requires 50+ hours of
#      wall-clock time. All 15 THP1 samples sharing one 72h SGE job is
#      structurally insufficient regardless of parallelism settings.
#      Fix: split into two phases:
#        Phase 1 — trax_map_sample (NEW, per-sample wildcard):
#          Each sample gets its own SGE job (1 slot, 72h budget). Runs the
#          bowtie2 + choosemappings.py pipeline directly, replicating what
#          processsamples.py does internally. Produces {sample}.bam in
#          trax/{cell_line}/ and a temp FASTQ in merged_bams/.
#          Per-sample BAMs and FASTQs are declared temp() — Snakemake
#          deletes them automatically after trax_quantify succeeds.
#        Phase 2 — trax_quantify (MODIFIED, per-cell-line as before):
#          Depends on all per-sample BAMs from Phase 1.
#          Runs processsamples.py --lazyremap, which detects the pre-built
#          BAMs in CWD and skips all bowtie2 mapping, running only the
#          counting, coverage, and R analysis steps. Observed < 2h for
#          A549; runtime budget set to 8h (480 min). Slots reduced to 8
#          (from 16) since parallel mapping is no longer needed here.
#      The rm -f *.bam pre-flight step is removed from trax_quantify —
#      those BAMs are now valid Snakemake-managed inputs, not stale files.
# =============================================================================

TRAX_CFG = config["trax"]
# NOTE: REF is defined in 00_reference_prep.smk (same included namespace).
#       Accessing config["references"] directly to avoid redefinition.


rule trax_map_sample:
    """
    FIX 10 Phase 1: per-sample mapping for TRAX.

    Runs the bowtie2 + choosemappings.py pipeline for ONE sample, replicating
    what processsamples.py does internally. Each sample gets its own 72h SGE
    job, so even the heaviest THP1 p8 samples (300M reads, ~55h bowtie2) fit
    within the wall-clock budget.

    Steps:
      1. Merge Pass 1 (functional) + Pass 2 (pre-tRNA) BAMs.
      2. Convert to FASTQ unconditionally (FIX 6: handles both single-end
         pass1 and paired-end pass2 reads; no -1/-2 flags).
      3. Delete merged BAM (pure intermediate; FASTQ is all TRAX needs).
      4. Write single-sample TRAX sample file.
      5. Run bowtie2 → choosemappings.py → samtools sort to produce
         {sample}.bam in trax/{cell_line}/ — exactly where
         processsamples.py --lazyremap will look for it.
      6. Index the BAM.

    Outputs are declared temp() so Snakemake deletes them automatically
    after trax_quantify (the downstream consumer) completes.
    """
    input:
        p1_bam  = f"{SCRATCH}/pass1_filters/{{sample}}/{{sample}}.functional.bam",
        p2_bam  = f"{SCRATCH}/pass2_pretRNA/{{sample}}/{{sample}}.pretRNA.bam",
        trax_db = config["references"]["trax_ref_dir"],
    output:
        # temp(): Snakemake deletes these once trax_quantify has consumed them.
        bam   = temp(f"{SCRATCH}/trax/{{cell_line}}/{{sample}}.bam"),
        bai   = temp(f"{SCRATCH}/trax/{{cell_line}}/{{sample}}.bam.bai"),
        fastq = temp(f"{SCRATCH}/trax/{{cell_line}}/merged_bams/{{sample}}.fastq"),
    params:
        outdir   = f"{SCRATCH}/trax/{{cell_line}}",
        trax_dir = config["trax"]["script_dir"],
    log:
        f"{SCRATCH}/logs/07_trax_map/{{cell_line}}/{{sample}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/07_trax_map/{{cell_line}}/{{sample}}.tsv",
    threads: 1
    resources:
        runtime   = config["resources"]["trax_map_sample"]["runtime_min"],
        sge_extra = sge_extra("trax_map_sample"),
    conda:
        "../../envs/trax_env.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/merged_bams
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        SAMPLE="{wildcards.sample}"
        CELL_LINE="{wildcards.cell_line}"
        DB="{input.trax_db}/hsapi38"
        MERGED="{params.outdir}/merged_bams/$SAMPLE.merged.bam"
        READS="{output.fastq}"
        OUTDIR="{params.outdir}"
        TRAX="{params.trax_dir}"

        echo "[$(date)] trax_map_sample: $SAMPLE (cell line: $CELL_LINE)"

        # --- Step 1: Merge pass1 + pass2 BAMs ---
        echo "[$(date)] Merging BAMs..."
        samtools merge -f -@ 1 "$MERGED" "{input.p1_bam}" "{input.p2_bam}"
        samtools index "$MERGED"

        # --- Step 2: Convert to FASTQ ---
        # FIX 6: plain samtools fastq — all reads output unconditionally.
        # Pass1 reads are single-end (no pairing flags); pass2 are paired-end.
        # bowtie2 maps with -U so pair order is irrelevant.
        echo "[$(date)] Converting to FASTQ..."
        samtools fastq -@ 1 "$MERGED" > "$READS"

        # --- Step 3: Clean up merged BAM (FASTQ is all TRAX needs) ---
        rm -f "$MERGED" "$MERGED.bai"

        # --- Step 4: Write single-sample TRAX sample file ---
        COND_CHAR=$(echo "$SAMPLE" | cut -d_ -f2 | cut -c1)
        COND="control"
        [ "$COND_CHAR" = "p" ] && COND="polyi"
        SAMPLE_FILE="$OUTDIR/merged_bams/$SAMPLE.sample.txt"
        printf "%s\t%s\t%s\n" "$SAMPLE" "$COND" "$READS" > "$SAMPLE_FILE"

        # --- Step 5: Map with bowtie2 + choosemappings.py ---
        # Replicates processsamples.py's internal per-sample mapping loop.
        # bowtie2 is hardcoded -p 1 by TRAX (no per-sample parallelism lever).
        # choosemappings.py routes reads to the sample BAM and filters
        # non-tRNA alignments below --minnontrnasize.
        # Use node-local scratch for sort temp if SGE set $TMPDIR (fast local
        # disk on Eddie compute nodes), otherwise fall back to outdir/tmp.
        # Written as explicit if/else rather than bash brace-expansion syntax
        # to avoid Snakemake's shell formatter misinterpreting dollar-braces.
        SORT_TMPDIR="$OUTDIR/tmp"
        if [ -n "$TMPDIR" ]; then SORT_TMPDIR="$TMPDIR"; fi
        mkdir -p "$SORT_TMPDIR"

        echo "[$(date)] Running bowtie2 + choosemappings.py..."
        cd "$OUTDIR"
        bowtie2 \
            -x "$DB-tRNAgenome" \
            -k 100 \
            --very-sensitive \
            --ignore-quals \
            --np 5 \
            --reorder \
            -p 1 \
            -U "$READS" | \
        python "$TRAX/choosemappings.py" \
            "$DB-trnatable.txt" \
            --progname=TRAX \
            --fqname="$READS" \
            --expname="$SAMPLE_FILE" \
            --minnontrnasize=20 | \
        samtools sort \
            -T "$SORT_TMPDIR/$SAMPLE-temp" \
            - \
            -o "$OUTDIR/$SAMPLE.bam"

        # --- Step 6: Index the BAM ---
        samtools index "$OUTDIR/$SAMPLE.bam"

        # --- Clean up single-sample file (not needed downstream) ---
        rm -f "$SAMPLE_FILE"

        echo "[$(date)] Mapping complete: $OUTDIR/$SAMPLE.bam"
        """


rule trax_quantify:
    """
    FIX 10 Phase 2: TRAX tRF quantification — analysis only (--lazyremap).

    All per-sample BAMs are pre-built by trax_map_sample. This rule runs
    processsamples.py --lazyremap, which detects the existing BAMs in the
    working directory and skips all bowtie2 mapping, running only the
    counting, coverage, and R analysis steps. Observed wall-clock < 2h for
    A549 (vs 23h total including mapping).

    The rm -f *.bam pre-flight step from the original rule is removed:
    those BAMs are now valid Snakemake-managed inputs (temp() outputs of
    trax_map_sample), not stale files to be discarded.

    tRAX scripts are called via:  python {params.trax_dir}/processsamples.py
    because tRAX is a GitHub repo, not a conda/pip package.  Its dependencies
    (samtools, pysam, bowtie2) are provided by trax_env.yaml.
    """
    input:
        # Per-sample BAMs and FASTQs from Phase 1 (trax_map_sample).
        # Declaring them here ensures trax_quantify waits for all per-sample
        # mapping to complete before running the analysis phase.
        # These are temp() outputs: Snakemake deletes them after this rule.
        per_sample_bams = lambda wildcards: expand(
            f"{SCRATCH}/trax/{{cell_line}}/{{sample}}.bam",
            cell_line=wildcards.cell_line,
            sample=samples_for(wildcards.cell_line)
        ),
        per_sample_fastqs = lambda wildcards: expand(
            f"{SCRATCH}/trax/{{cell_line}}/merged_bams/{{sample}}.fastq",
            cell_line=wildcards.cell_line,
            sample=samples_for(wildcards.cell_line)
        ),
        trax_db     = config["references"]["trax_ref_dir"],
        ensembl_gtf = config["references"]["ensembl_gtf"],
    output:
        # FIX: corrected filenames to match actual tRAX output (verified from
        # ls trax/A549/A549/ after a completed run). The originally declared
        # names (-readcounts.txt, -trnacounts.txt, -anticodoncounts.txt) do not
        # exist; tRAX produces these instead:
        readcounts      = f"{SCRATCH}/trax/{{cell_line}}/{{cell_line}}/{{cell_line}}-normalizedreadcounts.txt",
        trnacounts      = f"{SCRATCH}/trax/{{cell_line}}/{{cell_line}}/{{cell_line}}-trnamapinfo.txt",
        anticodoncounts = f"{SCRATCH}/trax/{{cell_line}}/{{cell_line}}/{{cell_line}}-aminocounts.txt",
        sample_file     = f"{SCRATCH}/trax/{{cell_line}}/trax_samples.txt",
    params:
        outdir     = f"{SCRATCH}/trax/{{cell_line}}",
        min_unique = TRAX_CFG["min_unique_cov"],
        # FIX 10: 8 cores for analysis phase (counting + R); the 16-core
        # setting was needed for concurrent mapping waves, which are now
        # handled per-sample by trax_map_sample.
        threads    = config["resources"]["trax_quantify"]["slots"],
        trax_dir   = config["trax"]["script_dir"],
        samples    = lambda wildcards: " ".join(samples_for(wildcards.cell_line)),
    log:
        f"{SCRATCH}/logs/07_trax/{{cell_line}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/07_trax/{{cell_line}}.tsv",
    threads: lambda wildcards: config["resources"]["trax_quantify"]["slots"]
    resources:
        runtime   = config["resources"]["trax_quantify"]["runtime_min"],
        sge_extra = sge_extra("trax_quantify"),
    conda:
        "../../envs/trax_env.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}/merged_bams
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        echo "[$(date)] TRAX analysis phase (--lazyremap) for {wildcards.cell_line}"

        # Recreate the full sample file pointing to FASTQs from trax_map_sample.
        # --lazyremap reads this file to know what samples to expect, then
        # checks for {{sample}}.bam in CWD and skips bowtie2 if found.
        > {output.sample_file}
        for SAMPLE in {params.samples}; do
            READS={params.outdir}/merged_bams/$SAMPLE.fastq
            COND_CHAR=$(echo "$SAMPLE" | cut -d_ -f2 | cut -c1)
            COND="control"
            [ "$COND_CHAR" = "p" ] && COND="polyi"
            printf "%s\t%s\t%s\n" "$SAMPLE" "$COND" "$READS" >> {output.sample_file}
        done

        echo "[$(date)] Sample file written:"
        cat {output.sample_file}

        # Run TRAX analysis only — all per-sample BAMs already exist in CWD.
        # --lazyremap: TRAX checks for {{sample}}.bam in CWD, skips mapping.
        # NOTE: rm -f *.bam pre-flight removed (FIX 10) — those BAMs are now
        # valid Snakemake-managed inputs from trax_map_sample, not stale files.
        echo "[$(date)] Running TRAX processsamples.py --lazyremap..."
        cd {params.outdir}
        python {params.trax_dir}/processsamples.py \
            --experimentname {wildcards.cell_line} \
            --databasename   {input.trax_db}/hsapi38 \
            --samplefile     {output.sample_file} \
            --mincoverage    {params.min_unique} \
            --cores          {params.threads} \
            --ensembl        {input.ensembl_gtf} \
            --lazyremap

        echo "[$(date)] TRAX complete for {wildcards.cell_line}."

        # FIX 7B: delete large merged_bams/ intermediates.
        # Per-sample BAMs (*.bam in this dir) are temp() outputs of
        # trax_map_sample — Snakemake deletes them automatically after this
        # rule succeeds; no manual rm needed for those.
        echo "[$(date)] Cleaning up merged_bams/ intermediates..."
        rm -f {params.outdir}/merged_bams/*.merged.bam
        rm -f {params.outdir}/merged_bams/*.merged.bam.bai
        rm -f {params.outdir}/merged_bams/*.fastq
        echo "[$(date)] Cleanup complete."
        """
