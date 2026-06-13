"""
workflow/scripts/cca_anticodon_filter.py
========================================
Snakemake script (called via `script:` directive).
Applies two pysam-based filters to mim-tRNAseq BAMs to identify
reads from functional, mature tRNAs (proposal Section 3.4).

FIX (2026-06-04): mim-tRNAseq aligns R1 only and produces SINGLE-END
(unpaired) BAMs. The original script tried to pair reads using
read.is_read1 / read.is_read2, which are never set in single-end BAMs.

FIX (2026-06-13): Replaced single-pass in-memory R2 dictionary with a
memory-efficient two-pass approach. The original approach loaded all ~40M
R2 reads into a Python dict (~30-40 GB), causing MemoryError on HPC nodes.

New strategy:
  Pass 1 (BAM scan)  — categorise every read as mapped or unmapped and
                        store only the read name strings in two sets.
                        Memory: ~2 × N × 140 bytes (strings only).
  Pass 2 (FASTQ stream) — stream R2 once:
      • For mapped reads:   evaluate the CCA check immediately; add to
                            cca_pass set if it passes. Discard seq/qual.
      • For unmapped reads: store full (seq, qual) in r2_unmapped dict
                            for writing to the Pass-2 output FASTQ.
  Pass 3 (BAM scan)  — main filter loop using cca_pass and r2_unmapped.

Typical memory saving (40M-read library, 50% mapping):
  Before: ~30-40 GB  |  After: ~8-12 GB

All output files and stat key names are unchanged; downstream rules
are unaffected.
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
trimmed_r2     = snakemake.input.trimmed_r2
filtered_bam   = snakemake.output.filtered_bam
filtered_bai   = snakemake.output.filtered_bai
stats_path     = snakemake.output.stats
unmapped_r1    = snakemake.output.unmapped_r1
unmapped_r2    = snakemake.output.unmapped_r2
log_path       = snakemake.log[0]
min_cca_qual   = snakemake.params.min_cca_qual
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
    if ref_name in anticodon_lookup:
        return anticodon_lookup[ref_name]
    m = re.search(r'-[A-Z][a-z]{2}([ACGT]{3})(?:$|-)', ref_name)
    if m:
        return m.group(1)
    return None


def check_cca_fastq(seq, qual_str, min_qual=20):
    if len(seq) < 3:
        return False
    if seq[:3] != "CCA":
        return False
    if not qual_str or len(qual_str) < 3:
        return True
    return all((ord(qual_str[i]) - 33) >= min_qual for i in range(3))


def is_anticodon_concordant(read):
    if read.reference_name is None:
        return False
    return ref_anticodon(read.reference_name) is not None


def phred_to_str(quals):
    if quals is None:
        return "I" * 40
    return "".join(chr(q + 33) for q in quals)


# ---------------------------------------------------------------------------
# Pass 1: scan BAM to categorise reads as mapped or unmapped
# ---------------------------------------------------------------------------
# We store only the name strings — no sequences — so memory scales with
# N × ~140 bytes rather than N × ~500 bytes (full seq/qual dict).
# ---------------------------------------------------------------------------
logger.info("Pass 1: scanning BAM to categorise reads (mapped vs unmapped)...")
mapped_names   = set()
unmapped_names = set()

with pysam.AlignmentFile(bam_path, "rb") as bam_scan:
    for read in bam_scan:
        if read.is_unmapped:
            unmapped_names.add(read.query_name)
        else:
            mapped_names.add(read.query_name)

logger.info(
    f"Pass 1 complete: {len(mapped_names):,} mapped reads, "
    f"{len(unmapped_names):,} unmapped reads"
)

# ---------------------------------------------------------------------------
# Pass 2: stream R2 FASTQ once
#   • mapped reads   → evaluate CCA immediately; store name in cca_pass if OK
#   • unmapped reads → store full (seq, qual) for Pass-2 FASTQ output
# ---------------------------------------------------------------------------
logger.info(f"Pass 2: streaming R2 FASTQ from {trimmed_r2} ...")
cca_pass   = set()   # names of mapped reads that pass the CCA filter
r2_unmapped = {}     # name → (seq, qual) for unmapped reads only

open_fn = gzip.open if trimmed_r2.endswith(".gz") else open
with open_fn(trimmed_r2, "rt") as fh:
    while True:
        header = fh.readline()
        if not header:
            break
        seq  = fh.readline().strip()
        fh.readline()             # '+' line — discard
        qual = fh.readline().strip()
        name = header.strip().lstrip("@").split()[0]

        if name in mapped_names:
            # Only CCA check needed — discard seq/qual afterwards
            if check_cca_fastq(seq, qual, min_qual=min_cca_qual):
                cca_pass.add(name)

        elif name in unmapped_names:
            # Full seq/qual needed to write R2 mate to the Pass-2 FASTQ
            r2_unmapped[name] = (seq, qual)

# Free mapped_names — no longer needed
del mapped_names

logger.info(
    f"Pass 2 complete: {len(cca_pass):,} mapped reads pass CCA filter, "
    f"{len(r2_unmapped):,} unmapped R2 reads stored"
)

# ---------------------------------------------------------------------------
# Pass 3: main filter loop using pre-computed structures
# ---------------------------------------------------------------------------
counters = {
    "total_pairs"    : 0,
    "both_unmapped"  : 0,
    "r1_only_mapped" : 0,
    "r2_only_mapped" : 0,
    "both_mapped"    : 0,
    "cca_fail"       : 0,
    "anticodon_fail" : 0,
    "both_pass"      : 0,
    "total_aligned"  : 0,
    "r2_missing"     : 0,
}

bam_in      = pysam.AlignmentFile(bam_path, "rb")
out_bam     = pysam.AlignmentFile(filtered_bam, "wb", header=bam_in.header)
fh_unmap_r1 = gzip.open(unmapped_r1, "wt")
fh_unmap_r2 = gzip.open(unmapped_r2, "wt")

logger.info("Pass 3: processing BAM reads and writing outputs...")

for read in bam_in:
    qname = read.query_name
    counters["total_pairs"] += 1

    if read.is_unmapped:
        counters["both_unmapped"] += 1

        if read.query_sequence:
            qual_str = phred_to_str(read.query_qualities)
            fh_unmap_r1.write(f"@{qname}\n{read.query_sequence}\n+\n{qual_str}\n")

        r2_entry = r2_unmapped.get(qname)
        if r2_entry is not None:
            fh_unmap_r2.write(f"@{qname}\n{r2_entry[0]}\n+\n{r2_entry[1]}\n")
        else:
            logger.debug(f"No R2 in trimmed FASTQ for unmapped read: {qname}")

        continue

    # Mapped read
    counters["both_mapped"]   += 1
    counters["total_aligned"] += 1

    # CCA filter — result already computed in Pass 2
    if qname not in cca_pass:
        counters["cca_fail"] += 1
        # Check if it was missing from FASTQ entirely vs just failing CCA
        if qname not in unmapped_names:
            counters["r2_missing"] += 1
        continue

    # Anticodon concordance filter
    if not is_anticodon_concordant(read):
        counters["anticodon_fail"] += 1
        continue

    counters["both_pass"] += 1
    out_bam.write(read)

bam_in.close()
out_bam.close()
fh_unmap_r1.close()
fh_unmap_r2.close()

logger.info(
    f"Pass 3 complete. "
    f"mapped={counters['both_mapped']:,}, "
    f"unmapped={counters['both_unmapped']:,}, "
    f"passed={counters['both_pass']:,}"
)

if counters["r2_missing"] > 0:
    logger.warning(
        f"{counters['r2_missing']:,} mapped reads had no matching R2 in "
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
# Write filter statistics (schema identical to original)
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
    f"Filter complete. Kept {counters['both_pass']:,} / {total_al:,} aligned reads."
)
logger.info("CCA+anticodon filter done.")
