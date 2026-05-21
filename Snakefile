# =============================================================================
# tRNA-seq Analysis Pipeline
# Project: Characterising tRNA library dynamics in the antiviral response
# Author:  Darren Lim
# =============================================================================
#
# Pipeline stages
# ---------------
#   00  Reference preparation  (pre-tRNA FASTA, Bowtie2 index, anticodon map,
#                                TRAX database, miRBase filtering)
#   01  Adapter trimming        (Trim Galore; on-the-fly A549 lane merging)
#   02  Post-trim QC            (FastQC + MultiQC per cell line)
#   03  Pass 1 alignment        (mim-tRNAseq → mature tRNA, per cell line)
#   04  Pass 1 filters          (pysam: CCA filter + anticodon concordance)
#   05  Pass 2 alignment        (Bowtie2 → pre-tRNA; unmapped from Pass 1)
#   06  Pass 3 screening        (miRDeep2 → miRNA; unmapped from Pass 2)
#   07  tRF quantification      (TRAX; Pass 1 + Pass 2 BAMs)
#   08  DESeq2 input prep       (count matrices + coldata per cell line)
#
# Usage
# -----
#   # Dry run (check DAG, no execution)
#   snakemake -n --cores 128
#
#   # Full run — activate the single conda env first, then:
#   conda activate snakemake_pipeline
#   snakemake --cores 128 --rerun-incomplete \
#             --keep-going --latency-wait 60
#
#   # Run in background (recommended on SSH server)
#   conda activate snakemake_pipeline
#   nohup snakemake --cores 128 --rerun-incomplete \
#         --keep-going --latency-wait 60 \
#         > logs/snakemake_main.log 2>&1 &
#   echo "Snakemake PID: $!"
#
# =============================================================================

import pandas as pd
import os

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
configfile: "config/config.yaml"

SCRATCH = config["scratch"]

# ---------------------------------------------------------------------------
# Load sample manifest
# ---------------------------------------------------------------------------
manifest = (
    pd.read_csv(config["manifest"], sep="\t", index_col="sample_id")
)

SAMPLES    = manifest.index.tolist()                          # all 30 samples
CELL_LINES = sorted(manifest["cell_line"].unique().tolist())  # ["A549", "THP1"]


def samples_for(cell_line):
    """Return list of sample_ids belonging to a given cell line."""
    return manifest[manifest["cell_line"] == cell_line].index.tolist()


def r_files(sample, read):
    """
    Return list of raw FASTQ paths for a sample and read (R1 or R2).
    Multi-lane A549 samples have 2 entries; THP1 has 1.
    """
    col  = "R1_files" if read == "R1" else "R2_files"
    return manifest.loc[sample, col].split(";")


def is_gzipped(sample):
    """True if the first R1 file for this sample is gzip-compressed."""
    first = r_files(sample, "R1")[0]
    return first.endswith(".gz")


def cat_cmd(sample, read):
    """
    Build a bash process-substitution expression that streams all lane files
    for this sample/read as a single uncompressed FASTQ stream.
    Used by Trim Galore to avoid pre-writing merged lane files to disk.

    A549 (.fastq, uncompressed):   <(cat  lane1.fastq lane2.fastq)
    THP1 (.fq.gz, gzip):           <(zcat file.fq.gz)
    """
    files   = r_files(sample, read)
    decat   = "zcat" if is_gzipped(sample) else "cat"
    joined  = " ".join(files)
    return f"<({decat} {joined})"


# ---------------------------------------------------------------------------
# Include rule modules (order is for readability only; Snakemake resolves DAG)
# ---------------------------------------------------------------------------
include: "workflow/rules/00_reference_prep.smk"
include: "workflow/rules/01_trim.smk"
include: "workflow/rules/02_post_trim_qc.smk"
include: "workflow/rules/03_pass1_mimtrnaseq.smk"
include: "workflow/rules/04_pass1_filters.smk"
include: "workflow/rules/05_pass2_pretRNA.smk"
include: "workflow/rules/06_pass3_mirna.smk"
include: "workflow/rules/07_trax.smk"
include: "workflow/rules/08_deseq2_prep.smk"

# ---------------------------------------------------------------------------
# Target rule  —  requesting all final outputs drives the entire DAG
# ---------------------------------------------------------------------------
rule all:
    input:
        # ── Post-trim QC ──────────────────────────────────────────────────
        expand(
            "{scratch}/qc/multiqc_post_trim/{cell_line}_post_trim_multiqc.html",
            scratch=SCRATCH,
            cell_line=CELL_LINES,
        ),
        # ── Pass 1: mim-tRNAseq isodecoder count tables ──────────────────
        expand(
            "{scratch}/pass1_mimtrnaseq/{cell_line}/counts/Isodecoder_counts.txt",
            scratch=SCRATCH,
            cell_line=CELL_LINES,
        ),
        # ── Pass 1: per-sample CCA+anticodon filtered BAMs ───────────────
        expand(
            "{scratch}/pass1_filters/{sample}/{sample}.functional.bam",
            scratch=SCRATCH,
            sample=SAMPLES,
        ),
        # ── Pass 2: per-sample pre-tRNA BAMs ─────────────────────────────
        expand(
            "{scratch}/pass2_pretRNA/{sample}/{sample}.pretRNA.bam",
            scratch=SCRATCH,
            sample=SAMPLES,
        ),
        # ── Pass 3: per-sample miRNA read counts ─────────────────────────
        expand(
            "{scratch}/pass3_mirna/{sample}/{sample}_miRNA_counts.txt",
            scratch=SCRATCH,
            sample=SAMPLES,
        ),
        # ── TRAX: per-cell-line tRF count tables ─────────────────────────
        expand(
            "{scratch}/trax/{cell_line}/counts/tRF_counts.txt",
            scratch=SCRATCH,
            cell_line=CELL_LINES,
        ),
        # ── DESeq2 inputs ─────────────────────────────────────────────────
        expand(
            "{scratch}/deseq2_input/{cell_line}/isodecoder_counts_matrix.tsv",
            scratch=SCRATCH,
            cell_line=CELL_LINES,
        ),
        expand(
            "{scratch}/deseq2_input/{cell_line}/isoacceptor_counts_matrix.tsv",
            scratch=SCRATCH,
            cell_line=CELL_LINES,
        ),
        expand(
            "{scratch}/deseq2_input/{cell_line}/coldata.tsv",
            scratch=SCRATCH,
            cell_line=CELL_LINES,
        ),
        # ── QC summary: read assignment proportions ───────────────────────
        expand(
            "{scratch}/qc/read_assignment/{cell_line}_read_assignment_summary.tsv",
            scratch=SCRATCH,
            cell_line=CELL_LINES,
        ),
