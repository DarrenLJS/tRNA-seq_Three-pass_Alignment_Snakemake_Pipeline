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
#   counts/Isoacceptor_counts.txt  → collapsed isoacceptor analysis (anticodon level)
#   mismatch/                      → wobble position 34 modification inference
#   align/{sample}.bam             → per-sample BAMs for downstream filters
#
# NOTE on mimseq output filenames (verified from actual run):
#   mimseq writes  counts/Isodecoder_counts_raw.txt   (NOT Isodecoder_counts.txt)
#                  counts/Anticodon_counts_raw.txt     (NOT Isoacceptor_counts.txt)
#   Also produced: counts/{cyto,mito}_Isodecoder_counts_DESEqNormalized.csv
#                  counts/{cyto,mito}_Anticodon_counts_DESEqNormalized.csv
#   We copy the raw files to the Snakemake-contracted names below.
#   The Anticodon → Isoacceptor renaming is intentional: anticodon groupings
#   are the isoacceptor level of analysis used throughout this pipeline.
#
# Conda environments
# ------------------
#   rule mimtrnaseq         → ../../envs/mimseq.yaml      (python=3.7, mimseq)
#   rule link_mimtrnaseq_bam → ../../envs/environment.yaml (samtools only)
# =============================================================================

MIM = config["mimtrnaseq"]
# NOTE: REF is already defined in 00_reference_prep.smk (same namespace).
#       Accessing config["references"] directly below to avoid redefinition.


def mim_input_r1(cell_line):
    """Return list of trimmed R1 FASTQ paths for all samples of a cell line."""
    return [
        f"{SCRATCH}/trimmed/{s}/{s}_val_1.fq.gz"
        for s in samples_for(cell_line)
    ]


rule mimtrnaseq:
    """
    Run mim-tRNAseq on all samples for one cell line.

    The current mimseq CLI (verified via --help) takes a tab-separated
    sampledata file as its positional argument:
        col 1: full path to FASTQ (R1 only for paired-end tRNA-seq)
        col 2: condition/group name

    Conditions are extracted from the trimmed FASTQ filenames:
        *_c[0-9]* -> "c"  (control)
        *_p[0-9]* -> "p"  (perturbation/treatment)

    IMPORTANT: mimseq requires --out-dir to be a non-existing directory.
    Snakemake pre-creates {params.outdir} for the sentinel/count outputs,
    so we point mimseq at a fresh _run/ subdirectory and then copy the
    required outputs up to the Snakemake-expected paths afterwards.

    Per-sample BAMs land at {outdir}/align/{sample}_val_1.fq.gz.unpaired_uniq.bam
    (mimseq preserves the full input filename and appends the GSNAP category suffix),
    which is where link_mimtrnaseq_bam expects to find them.
    """
    input:
        r1_files = lambda wildcards: mim_input_r1(wildcards.cell_line),
    output:
        iso_counts    = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isodecoder_counts.txt",
        isoa_counts   = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isoacceptor_counts.txt",
        align_done    = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/.align_done",
        mismatch_done = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/.mismatch_done",
    params:
        outdir       = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}",
        # mimseq writes here; must not pre-exist (Snakemake creates outdir above)
        mimseq_dir   = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/_run",
        species      = config["references"]["mimtrnaseq_species"],  # 'Hsap' for hg38
        threads      = MIM["threads"],
        min_cov      = MIM["min_cov"],
        max_multi    = MIM["max_multi"],
        cluster_id   = MIM["cluster_id"],
        control_cond = "c",                                          # control condition label
        name         = lambda wildcards: f"{wildcards.cell_line}_tRNAseq",
        r1_csv       = lambda wildcards, input: ",".join(input.r1_files),
    log:
        f"{SCRATCH}/logs/03_mimtrnaseq/{{cell_line}}.log",
    benchmark:
        f"{SCRATCH}/benchmarks/03_mimtrnaseq/{{cell_line}}.tsv",
    threads: lambda wildcards: MIM["threads"]
    resources:
        # FIX: was sge_pe="sharedmem" + a hardcoded "-V -l h_vmem=8000M"
        # literal — same broken pattern already fixed in trim_galore,
        # bowtie2_pretRNA, and trax_quantify (sge_pe collapses to 1 SGE
        # slot under the EDDIE profile's --cores 1, regardless of
        # threads:). mim-tRNAseq's GSNAP alignment genuinely uses
        # --threads {params.threads} (8) in parallel, so this was a real
        # latent OOM risk that simply hadn't been hit yet. sge_extra()
        # hardcodes "-pe sharedmem N" so the plugin can't collapse it;
        # slots/vmem/runtime now live in config["resources"]["mimtrnaseq"].
        runtime   = config["resources"]["mimtrnaseq"]["runtime_min"],
        sge_extra = sge_extra("mimtrnaseq"),
    conda:
        "../../envs/mimseq.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Starting mim-tRNAseq for cell line: {wildcards.cell_line}" > {log} 2>&1

        # ── Build sample data sheet ──────────────────────────────────────
        # mimseq takes a tab-separated file: col1=fastq_path, col2=condition.
        # Condition is extracted from the filename (_c[digit] or _p[digit]).
        SAMPLESHEET="$(dirname {log})/{wildcards.cell_line}_samplesheet.txt"
        > "$SAMPLESHEET"
        for f in $(echo "{params.r1_csv}" | tr ',' '\n'); do
            COND=$(basename "$f" | sed 's/.*_\([cp]\)[0-9].*/\1/')
            printf "%s\t%s\n" "$f" "$COND" >> "$SAMPLESHEET"
        done
        echo "[$(date)] Sample sheet:" >> {log}
        cat "$SAMPLESHEET" >> {log}

        # ── Run mim-tRNAseq ──────────────────────────────────────────────
        # Use a fresh _run/ subdir; mimseq requires --out-dir to not exist.
        rm -rf "{params.mimseq_dir}"

        # Verify real usearch is on PATH.
        # mimseq calls `usearch -cluster_fast` then `usearch -sortbysize`.
        # vsearch does NOT add ;size= annotations during clustering, so the
        # subsequent sortbysize call fails with exit code 1.  A vsearch symlink
        # is therefore insufficient — usearch v11 must be installed in the env:
        #   wget https://drive5.com/downloads/usearch11.0.667_i86linux32.gz
        #   gunzip usearch11.0.667_i86linux32.gz
        #   chmod +x usearch11.0.667_i86linux32
        #   cp usearch11.0.667_i86linux32 $CONDA_PREFIX/bin/usearch
        if ! command -v usearch &>/dev/null; then
            echo "[ERROR] usearch not found on PATH. Install usearch v11 in the mimseq conda env." >> {log}
            echo "[ERROR] See comment in 03_pass1_mimtrnaseq.smk for installation instructions." >> {log}
            exit 1
        fi
        echo "[$(date)] Using usearch: $(which usearch)" >> {log}

        mimseq \
            --species           {params.species} \
            --cluster-id        {params.cluster_id} \
            --threads           {params.threads} \
            --min-cov           {params.min_cov} \
            --max-multi         {params.max_multi} \
            --control-condition {params.control_cond} \
            --local-mod \
            -n                  {params.name} \
            --out-dir           "{params.mimseq_dir}" \
            "$SAMPLESHEET" \
            >> {log} 2>&1

        # ── Move outputs to Snakemake-expected paths ─────────────────────
        # mimseq actual output names (verified from run):
        #   Isodecoder_counts_raw.txt  → Isodecoder_counts.txt  (isodecoder level)
        #   Anticodon_counts_raw.txt   → Isoacceptor_counts.txt (isoacceptor/anticodon level)
        # The Anticodon→Isoacceptor rename is intentional: anticodon groupings
        # represent the isoacceptor level used throughout this pipeline.
        echo "[$(date)] Copying count files..." >> {log}
        mkdir -p "$(dirname {output.iso_counts})"
        cp "{params.mimseq_dir}/counts/Isodecoder_counts_raw.txt" "{output.iso_counts}"
        cp "{params.mimseq_dir}/counts/Anticodon_counts_raw.txt"  "{output.isoa_counts}"

        # Move align/ to {params.outdir}/align/ for link_mimtrnaseq_bam
        echo "[$(date)] Moving align directory..." >> {log}
        rm -rf "{params.outdir}/align"
        mv "{params.mimseq_dir}/align" "{params.outdir}/align"

        # ── Verify and create sentinel files ─────────────────────────────
        N_BAMS=$(find "{params.outdir}/align/" -name "*.bam" | wc -l)
        echo "[$(date)] Found $N_BAMS BAM files in {params.outdir}/align/" >> {log}

        touch {output.align_done}
        touch {output.mismatch_done}
        echo "[$(date)] mim-tRNAseq complete for {wildcards.cell_line}." >> {log}
        """


# ---------------------------------------------------------------------------
# Helper rule: locate the per-sample BAM produced by mim-tRNAseq.
# mim-tRNAseq preserves the full input FASTQ filename and appends
# '.unpaired_uniq.bam' (GSNAP uniquely-mapping output category), e.g.:
#   A549_c2_1_val_1.fq.gz.unpaired_uniq.bam
# No unpaired_mult BAMs are present — mimseq handles multimapper filtering
# internally. We copy to our consistent per-sample naming scheme so downstream
# rules can use {sample}.bam paths without worrying about the mim naming logic.
# ---------------------------------------------------------------------------
rule link_mimtrnaseq_bam:
    """
    Create a consistently-named copy of each per-sample mim-tRNAseq BAM.
    mim-tRNAseq preserves the full input FASTQ filename and appends
    '.unpaired_uniq.bam' (e.g. A549_c2_1_val_1.fq.gz.unpaired_uniq.bam).
    We copy to {sample}.bam for all downstream rules.

    Uses the main environment (samtools only; mimseq not required here).
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
        # mim-tRNAseq preserves the full input FASTQ filename and appends
        # '.unpaired_uniq.bam' — confirmed from actual run output:
        #   A549_c2_1_val_1.fq.gz.unpaired_uniq.bam
        # Only unpaired_uniq BAMs are present; no unpaired_mult BAMs exist,
        # meaning mimseq handles multimapper filtering internally.
        src_bam = lambda wildcards: (
            f"{SCRATCH}/pass1_mimtrnaseq/"
            f"{manifest.loc[wildcards.sample,'cell_line']}/"
            f"align/{wildcards.sample}_val_1.fq.gz.unpaired_uniq.bam"
        ),
        outdir = f"{SCRATCH}/pass1_mimtrnaseq/{{sample}}",
    log:
        f"{SCRATCH}/logs/03_link_bam/{{sample}}.log",
    resources:
        runtime = 30,
    conda:
        "../../envs/environment.yaml"
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
