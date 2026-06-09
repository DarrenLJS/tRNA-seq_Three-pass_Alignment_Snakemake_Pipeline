"""
workflow/scripts/cca_anticodon_filter.py
========================================
Snakemake script (called via `script:` directive).
Applies two pysam-based filters to mim-tRNAseq BAMs to identify
reads from functional, mature tRNAs (proposal Section 3.4).

FIX (2026-06-04): mim-tRNAseq aligns R1 only and produces SINGLE-END
(unpaired) BAMs. The original script tried to pair reads using
read.is_read1 / read.is_read2, which are never set in single-end BAMs.
Every read therefore ended up as [None, None], every read was counted as
both_unmapped, and both unmapped output FASTQs were written empty —
causing all downstream bowtie2_pretRNA jobs to fail with "0 reads".

This rewrite iterates reads individually (no pairing by flag). The
original trimmed R2 FASTQ is loaded into a dict at startup so that:
  (a) the CCA check can be applied to R2 for mapped reads, and
  (b) R2 mates can be recovered for unmapped reads passed to Pass 2.

All output files and stat key names are unchanged; downstream rules
are unaffected.

Filter 1 — CCA filter
    Mature tRNAs carry a non-encoded 3'-CCA tail.
    R2 primes from the 3' end; the first three bases of R2 should be "CCA".
    R2 sequences are loaded from the original trimmed R2 FASTQ.

Filter 2 — Anticodon concordance filter
    The anticodon is extracted from the reference sequence name
    (e.g. chr1.tRNA1-AlaAGC -> anticodon AGC) and cross-checked against
    the anticodon map TSV.

Inputs (from snakemake.input):
    bam           : mim-tRNAseq aligned BAM (single-end, indexed)
    anticodon_map : TSV of locus -> anticodon (built by rule build_anticodon_map)
    trimmed_r2    : original trimmed R2 FASTQ (.fq.gz) for this sample

Outputs (from snakemake.output):
    filtered_bam  : BAM with reads passing BOTH filters (single-end)
    filtered_bai  : index for filtered_bam
    stats         : TSV of filter statistics (same schema as before)
    unmapped_r1   : R1 FASTQ for reads unmapped in Pass 1 (for Pass 2)
    unmapped_r2   : R2 FASTQ for reads unmapped in Pass 1 (for Pass 2)
"""

import pysam
import gzip
import logging
import os
import re

# ---------------------------------------------------------------------------
# Snakemake bindings
# ---------------------------------------------------------------------------
bam_path       = snakemake.input.bam
anticodon_map  = snakemake.input.anticodon_map
trimmed_r2     = snakemake.input.trimmed_r2    # FIX: original trimmed R2 FASTQ
filtered_bam   = snakemake.output.filtered_bam
filtered_bai   = snakemake.output.filtered_bai
stats_path     = snakemake.output.stats
unmapped_r1    = snakemake.output.unmapped_r1
unmapped_r2    = snakemake.output.unmapped_r2
log_path       = snakemake.log[0]
min_cca_qual   = snakemake.params.min_cca_qual   # default 20
sample_id      = snakemake.wildcards.sample

os.makedirs(snakemake.params.outdir, exist_ok=True)
os.makedirs(os.path.dirname(log_path), exist_ok=True)

logging.basicConfig(
    filename=log_path, level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger()
logger.info(f"Starting CCA+anticodon filter for sample: {sample_id}")
logger.info("Mode: single-end BAM (mim-tRNAseq R1-only output)")

# ---------------------------------------------------------------------------
# Load anticodon lookup table
# ---------------------------------------------------------------------------
anticodon_lookup = {}
with open(anticodon_map) as fh:
    for line in fh:
        if line.startswith("locus"):
            continue
        parts = line.strip().split("\t")
        if len(parts) == 2:
            anticodon_lookup[parts[0]] = parts[1]

logger.info(f"Loaded {len(anticodon_lookup)} locus->anticodon entries")


def ref_anticodon(ref_name):
    """
    Extract the expected anticodon from a tRNA reference sequence name.
    Two sources:
      1. anticodon_lookup dict (exact locus match)
      2. Parse from reference name: 'chr1.tRNA1-AlaAGC' -> 'AGC' (last 3 chars)
    Falls back to None if neither works.
    """
    if ref_name in anticodon_lookup:
        return anticodon_lookup[ref_name]
    m = re.search(r'-[A-Z][a-z]{2}([ACGT]{3})(?:$|-)', ref_name)
    if m:
        return m.group(1)
    return None


def check_cca_fastq(seq, qual_str, min_qual=20):
    """
    Return True if the first three bases of seq are 'CCA' and all three
    Phred+33-encoded base qualities meet min_qual.
    """
    if len(seq) < 3:
        return False
    if seq[:3] != "CCA":
        return False
    if not qual_str or len(qual_str) < 3:
        return True   # no quality info; accept on sequence alone
    return all((ord(qual_str[i]) - 33) >= min_qual for i in range(3))


def is_anticodon_concordant(read):
    """
    Return True if the reference name encodes a recognised anticodon.

    TODO: Implement position-level anticodon check using mim-tRNAseq
          reference coordinates once the alignment structure is confirmed.
    """
    if read.reference_name is None:
        return False
    return ref_anticodon(read.reference_name) is not None


def phred_to_str(quals):
    if quals is None:
        return "I" * 40   # dummy quality
    return "".join(chr(q + 33) for q in quals)


# ---------------------------------------------------------------------------
# Load R2 sequences from the trimmed FASTQ into memory
# ---------------------------------------------------------------------------
# Keys are read names (FASTQ header minus '@' and any trailing description).
# Values are (sequence_str, quality_str) tuples.
# Memory: ~12M reads x ~100 bytes = 1-2 GB; acceptable on HPC nodes.
logger.info(f"Loading R2 sequences from {trimmed_r2} ...")
r2_dict = {}
open_fn = gzip.open if trimmed_r2.endswith(".gz") else open
with open_fn(trimmed_r2, "rt") as fh:
    while True:
        header = fh.readline()
        if not header:
            break
        seq  = fh.readline().strip()
        fh.readline()             # '+' line — discard
        qual = fh.readline().strip()
        # Strip '@' and take only the first whitespace-delimited word so
        # that Illumina description fields (" 2:N:0:ATCACG") are ignored.
        name = header.strip().lstrip("@").split()[0]
        r2_dict[name] = (seq, qual)

logger.info(f"Loaded {len(r2_dict)} R2 reads into memory")

# ---------------------------------------------------------------------------
# Main filter loop — single-end reads
# ---------------------------------------------------------------------------
# Counter keys are kept identical to the original schema so that the
# build_read_assignment_summary rule and any QC reports are unaffected.
counters = {
    "total_pairs"    : 0,   # = total R1 reads processed (one per BAM record)
    "both_unmapped"  : 0,   # R1 reads that did not map in Pass 1
    "r1_only_mapped" : 0,   # unused in single-end mode; kept for schema compat
    "r2_only_mapped" : 0,   # unused in single-end mode; kept for schema compat
    "both_mapped"    : 0,   # R1 reads that mapped (before CCA/anticodon checks)
    "cca_fail"       : 0,
    "anticodon_fail" : 0,
    "both_pass"      : 0,
    "total_aligned"  : 0,
    "r2_missing"     : 0,   # mapped reads with no R2 in trimmed FASTQ (warning)
}

# FIX (memory): open unmapped FASTQ output handles BEFORE the BAM loop so
# reads are streamed directly to disk rather than accumulated in RAM lists.
# Previously unmapped_reads_r1/r2 were Python lists that held all unmapped
# reads in memory until after the loop, doubling peak RAM usage on top of
# the r2_dict that is already ~1–2 GB.
bam_in      = pysam.AlignmentFile(bam_path, "rb")
out_bam     = pysam.AlignmentFile(filtered_bam, "wb", header=bam_in.header)
fh_unmap_r1 = gzip.open(unmapped_r1, "wt")
fh_unmap_r2 = gzip.open(unmapped_r2, "wt")

logger.info("Processing BAM reads (single-end mode)...")

for read in bam_in:
    qname = read.query_name
    counters["total_pairs"] += 1

    # ------------------------------------------------------------------
    # Unmapped reads -> stream directly to Pass 2 FASTQs
    # ------------------------------------------------------------------
    if read.is_unmapped:
        counters["both_unmapped"] += 1

        # R1: sequence is preserved as-is in the BAM for unmapped reads
        if read.query_sequence:
            qual_str = phred_to_str(read.query_qualities)
            fh_unmap_r1.write(f"@{qname}\n{read.query_sequence}\n+\n{qual_str}\n")

        # R2: recover from the original trimmed FASTQ
        r2_entry = r2_dict.get(qname)
        if r2_entry is not None:
            fh_unmap_r2.write(f"@{qname}\n{r2_entry[0]}\n+\n{r2_entry[1]}\n")
        else:
            logger.debug(f"No R2 in trimmed FASTQ for unmapped read: {qname}")

        continue

    # ------------------------------------------------------------------
    # Mapped reads -> apply CCA + anticodon filters
    # ------------------------------------------------------------------
    counters["both_mapped"]   += 1
    counters["total_aligned"] += 1

    # CCA filter: look up R2 from the trimmed FASTQ
    r2_entry = r2_dict.get(qname)
    if r2_entry is None:
        counters["r2_missing"] += 1
        counters["cca_fail"]   += 1
        logger.debug(f"No R2 in trimmed FASTQ for mapped read: {qname}")
        continue

    r2_seq, r2_qual = r2_entry
    if not check_cca_fastq(r2_seq, r2_qual, min_qual=min_cca_qual):
        counters["cca_fail"] += 1
        continue

    # Anticodon concordance filter
    if not is_anticodon_concordant(read):
        counters["anticodon_fail"] += 1
        continue

    # Both filters passed -> write to output BAM
    counters["both_pass"] += 1
    out_bam.write(read)

bam_in.close()
out_bam.close()
fh_unmap_r1.close()
fh_unmap_r2.close()

logger.info(
    f"BAM iteration complete. "
    f"mapped={counters['both_mapped']}, "
    f"unmapped={counters['both_unmapped']}, "
    f"passed={counters['both_pass']}"
)
logger.info(
    f"Streamed {counters['both_unmapped']} unmapped R1 reads to {unmapped_r1}"
)
if counters["r2_missing"] > 0:
    logger.warning(
        f"{counters['r2_missing']} mapped reads had no matching R2 in "
        f"the trimmed FASTQ — verify that trimmed_r2 is for the correct sample."
    )

# ---------------------------------------------------------------------------
# Sort and index the filtered BAM
# ---------------------------------------------------------------------------
logger.info("Sorting and indexing filtered BAM...")
sorted_tmp = filtered_bam + ".sorted.tmp"
pysam.sort("-o", sorted_tmp, filtered_bam)
os.rename(sorted_tmp, filtered_bam)
pysam.index(filtered_bam)

# ---------------------------------------------------------------------------
# Write filter statistics  (schema identical to original)
# ---------------------------------------------------------------------------
total_al        = counters["total_aligned"]
functional_rate = counters["both_pass"] / total_al if total_al > 0 else 0.0
cca_rate        = (
    (total_al - counters["cca_fail"]) / total_al if total_al > 0 else 0.0
)

stats_rows = {
    "sample"         : sample_id,
    "total_pairs"    : counters["total_pairs"],
    "both_unmapped"  : counters["both_unmapped"],
    "both_mapped"    : counters["both_mapped"],
    "total_aligned"  : counters["total_aligned"],
    "cca_fail"       : counters["cca_fail"],
    "anticodon_fail" : counters["anticodon_fail"],
    "both_pass"      : counters["both_pass"],
    "cca_pass_rate"  : round(cca_rate, 4),
    "functional_rate": round(functional_rate, 4),
    "unmapped"       : counters["both_unmapped"],
}

with open(stats_path, "w") as sf:
    sf.write("metric\tvalue\n")
    for k, v in stats_rows.items():
        sf.write(f"{k}\t{v}\n")

# QC flag
thresh = 0.70
if functional_rate < thresh:
    logger.warning(
        f"SAMPLE FLAG: {sample_id} functional rate {functional_rate:.1%} "
        f"is below target {thresh:.0%}"
    )
else:
    logger.info(
        f"Functional rate: {functional_rate:.1%} (target >={thresh:.0%}) -- OK"
    )

logger.info(
    f"Filter complete. Kept {counters['both_pass']} / {total_al} aligned reads."
)
logger.info("CCA+anticodon filter done.")
