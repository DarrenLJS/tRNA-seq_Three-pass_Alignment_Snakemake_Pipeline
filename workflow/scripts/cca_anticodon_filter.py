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

FIX (2026-06-14) Bug 1 — CCA orientation:
    R2 is sequenced FROM the 3′ CCA tail of the tRNA in REVERSE COMPLEMENT
    orientation. CCA (5′→3′ on the tRNA) appears as TGG at the 5′ start of
    R2, not CCA. Confirmed by diagnostic: 34% of R2 reads start with TGG vs
    0.17% starting with CCA. Illumina R2 position 0 is also frequently
    miscalled as N, so NGG is also accepted. Quality check skips position 0.

FIX (2026-06-14) Bug 2 — Unmapped read recovery:
    mim-tRNAseq writes ONLY mapped reads to its BAM (100% mapping rate
    confirmed; 0 is_unmapped reads). The previous code found no unmapped
    reads in the BAM and wrote empty unmapped FASTQs, leaving ~55 M reads
    stranded with nothing feeding into Pass 2 (Bowtie2 pre-tRNA).
    Fix: stream trimmed_r1 and trimmed_r2 simultaneously; any read pair
    whose name is absent from mapped_names is written directly to disk
    without buffering in memory. trimmed_r1 is now a required input.

New three-pass strategy:
  Pass 1 (BAM scan)     — collect all mapped read names into mapped_names.
  Pass 2 (FASTQ stream) — stream R1 + R2 simultaneously:
      • name in mapped_names  → TGG/NGG CCA check on R2 → add to cca_pass
      • name not in mapped_names → write to unmapped_r1 / unmapped_r2
  Pass 3 (BAM scan)     — CCA + anticodon filter → write filtered_bam.

FIX (2026-07-02) — independent filter accounting:
    Pass 3 previously short-circuited: a CCA failure skipped the anticodon
    check entirely via `continue`, so `anticodon_fail` only ever meant
    "passed CCA, failed anticodon" and there was no way to tell "failed
    CCA only" apart from "failed both filters". Both checks now run on
    every read unconditionally, giving four independent buckets:
        both_pass            — passed CCA AND anticodon (written to output)
        fail_cca_only        — failed CCA, would have passed anticodon
        fail_anticodon_only  — passed CCA, failed anticodon
        fail_both            — failed both filters
    This also fixes a mislabeling bug where filter_stats.tsv's
    "anticodon_pass" column was silently set to bam_counters["both_pass"]
    (a duplicate of both_pass) instead of the true anticodon-pass count.
    filter_stats.tsv's schema is NOT identical to the original as a
    result — see the Pass 3 / stats_rows sections below for the new
    columns (fail_cca_only, fail_anticodon_only, fail_both) and the
    corrected anticodon_pass values.
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
trimmed_r1     = snakemake.input.trimmed_r1    # FIX (Bug 2): recover unmapped R1 reads absent from BAM
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


def check_cca_r2(seq, qual_str, min_qual=20):
    """
    FIX (Bug 1): R2 reads from the 3' CCA end of tRNA in REVERSE COMPLEMENT.
    CCA on the tRNA (5'→3') appears as TGG at the 5' start of R2.

    Two Illumina R2 artefacts handled:
      - Position 0 is frequently miscalled as N (low first-cycle confidence).
        Accept T or N at position 0.
      - Quality at position 0 is often '#' (Phred 2). Skip quality check there.
    """
    if len(seq) < 3:
        return False
    if seq[0] not in ("T", "N"):   # TGG or NGG both accepted
        return False
    if seq[1:3] != "GG":
        return False
    # Quality check on positions 1 and 2 only (skip unreliable pos 0)
    if qual_str and len(qual_str) >= 3:
        return all((ord(qual_str[i]) - 33) >= min_qual for i in range(1, 3))
    return True


def is_anticodon_concordant(read):
    if read.reference_name is None:
        return False
    return ref_anticodon(read.reference_name) is not None


def phred_to_str(quals):
    if quals is None:
        return "I" * 40
    return "".join(chr(q + 33) for q in quals)


# ---------------------------------------------------------------------------
# Pass 1: scan BAM to build mapped_names set
# ---------------------------------------------------------------------------
# FIX (Bug 2): mim-tRNAseq writes ONLY mapped reads to its BAM (100% mapping
# rate confirmed). There are no is_unmapped reads in this file. Unmapped reads
# are recovered by FASTQ subtraction in Pass 2.
# ---------------------------------------------------------------------------
logger.info("Pass 1: building mapped_names set from BAM (all reads are mapped)...")
mapped_names = set()

with pysam.AlignmentFile(bam_path, "rb") as bam_scan:
    for read in bam_scan:
        mapped_names.add(read.query_name)

logger.info(f"Pass 1 complete: {len(mapped_names):,} mapped read names")

# ---------------------------------------------------------------------------
# Pass 2: stream trimmed R1 + R2 FASTQs simultaneously
#
# FIX (Bug 1): CCA check now looks for TGG/NGG at R2 start (RC of CCA).
# FIX (Bug 2): unmapped reads are recovered here by FASTQ subtraction —
#   reads absent from mapped_names are written directly to disk without
#   buffering in memory (avoids storing ~55 M sequences as a Python dict).
# ---------------------------------------------------------------------------
logger.info(f"Pass 2: streaming R1+R2 FASTQs to check CCA and recover unmapped reads...")

cca_pass = set()   # read names whose R2 starts with TGG/NGG (= RC of CCA)

fastq_counters = {
    "total_input"     : 0,
    "mapped_seen"     : 0,
    "cca_pass"        : 0,
    "cca_fail"        : 0,
    "unmapped_written": 0,
    "name_mismatch"   : 0,
}

open_r1 = gzip.open if trimmed_r1.endswith(".gz") else open
open_r2 = gzip.open if trimmed_r2.endswith(".gz") else open

with (open_r1(trimmed_r1, "rt") as r1_fh,
      open_r2(trimmed_r2, "rt") as r2_fh,
      gzip.open(unmapped_r1, "wt") as out_r1,
      gzip.open(unmapped_r2, "wt") as out_r2):

    while True:
        # R1 record
        r1_hdr = r1_fh.readline()
        if not r1_hdr:
            break
        r1_seq  = r1_fh.readline().strip()
        r1_fh.readline()            # '+' line — discard
        r1_qual = r1_fh.readline().strip()
        r1_name = r1_hdr.strip().lstrip("@").split()[0]

        # R2 record
        r2_hdr  = r2_fh.readline()
        if not r2_hdr:
            logger.warning("R2 FASTQ ended before R1 — mismatched read counts")
            break
        r2_seq  = r2_fh.readline().strip()
        r2_fh.readline()            # '+' line — discard
        r2_qual = r2_fh.readline().strip()
        r2_name = r2_hdr.strip().lstrip("@").split()[0]

        fastq_counters["total_input"] += 1

        if r1_name != r2_name:
            fastq_counters["name_mismatch"] += 1
            logger.warning(
                f"R1/R2 name mismatch at pair {fastq_counters['total_input']:,}: "
                f"R1={r1_name}, R2={r2_name} — skipping pair"
            )
            continue

        if r1_name in mapped_names:
            fastq_counters["mapped_seen"] += 1
            if check_cca_r2(r2_seq, r2_qual, min_qual=min_cca_qual):
                cca_pass.add(r1_name)
                fastq_counters["cca_pass"] += 1
            else:
                fastq_counters["cca_fail"] += 1
        else:
            # Unmapped: write both mates for Pass 2 (Bowtie2 pre-tRNA)
            out_r1.write(f"@{r1_name}\n{r1_seq}\n+\n{r1_qual}\n")
            out_r2.write(f"@{r2_name}\n{r2_seq}\n+\n{r2_qual}\n")
            fastq_counters["unmapped_written"] += 1

# mapped_names no longer needed — free memory before Pass 3
del mapped_names

logger.info(
    f"Pass 2 complete: "
    f"total_input={fastq_counters['total_input']:,}, "
    f"mapped={fastq_counters['mapped_seen']:,}, "
    f"cca_pass={fastq_counters['cca_pass']:,} "
    f"({fastq_counters['cca_pass']/max(fastq_counters['mapped_seen'],1)*100:.1f}%), "
    f"unmapped_written={fastq_counters['unmapped_written']:,}"
)
if fastq_counters["name_mismatch"] > 0:
    logger.warning(
        f"{fastq_counters['name_mismatch']:,} R1/R2 name mismatches — "
        f"verify trimmed_r1 and trimmed_r2 are from the same sample."
    )

# ---------------------------------------------------------------------------
# Pass 3: scan BAM and apply CCA + anticodon filters → write filtered_bam
# ---------------------------------------------------------------------------
# Unmapped reads are already written to disk in Pass 2.
# This pass only handles mapped reads: CCA result from cca_pass, then
# anticodon concordance check.
# ---------------------------------------------------------------------------
bam_counters = {
    "total_aligned"      : 0,
    "cca_pass"           : 0,
    "cca_fail"           : 0,
    "anticodon_pass"     : 0,
    "anticodon_fail"     : 0,
    "fail_cca_only"      : 0,   # failed CCA, would have passed anticodon
    "fail_anticodon_only": 0,   # passed CCA, failed anticodon
    "fail_both"          : 0,   # failed both filters
    "both_pass"          : 0,   # passed both filters (written to output BAM)
}

bam_in  = pysam.AlignmentFile(bam_path, "rb")
out_bam = pysam.AlignmentFile(filtered_bam, "wb", header=bam_in.header)

logger.info("Pass 3: filtering BAM reads (CCA + anticodon, checked independently)...")

for read in bam_in:
    bam_counters["total_aligned"] += 1

    # FIX (read-count reporting): both checks now run unconditionally on
    # every read, rather than short-circuiting after a CCA failure. This
    # is required to distinguish "failed CCA only" from "failed both
    # filters" — the previous version could not tell these apart, since
    # a CCA failure skipped the anticodon check entirely.
    cca_ok       = read.query_name in cca_pass
    anticodon_ok = is_anticodon_concordant(read)

    if cca_ok:
        bam_counters["cca_pass"] += 1
    else:
        bam_counters["cca_fail"] += 1

    if anticodon_ok:
        bam_counters["anticodon_pass"] += 1
    else:
        bam_counters["anticodon_fail"] += 1

    if cca_ok and anticodon_ok:
        bam_counters["both_pass"] += 1
        out_bam.write(read)
    elif cca_ok and not anticodon_ok:
        bam_counters["fail_anticodon_only"] += 1
    elif not cca_ok and anticodon_ok:
        bam_counters["fail_cca_only"] += 1
    else:
        bam_counters["fail_both"] += 1

bam_in.close()
out_bam.close()

logger.info(
    f"Pass 3 complete: "
    f"total_aligned={bam_counters['total_aligned']:,}, "
    f"cca_pass={bam_counters['cca_pass']:,}, cca_fail={bam_counters['cca_fail']:,}, "
    f"anticodon_pass={bam_counters['anticodon_pass']:,}, anticodon_fail={bam_counters['anticodon_fail']:,}, "
    f"fail_cca_only={bam_counters['fail_cca_only']:,}, "
    f"fail_anticodon_only={bam_counters['fail_anticodon_only']:,}, "
    f"fail_both={bam_counters['fail_both']:,}, "
    f"both_pass={bam_counters['both_pass']:,}"
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
# Write filter statistics (schema CHANGED from original — see FIX 2026-07-02
# in the module docstring: new fail_cca_only/fail_anticodon_only/fail_both
# columns, and anticodon_pass is now the true anticodon-pass count rather
# than a duplicate of both_pass)
# ---------------------------------------------------------------------------
total_al        = bam_counters["total_aligned"]
functional_rate = bam_counters["both_pass"] / total_al if total_al > 0 else 0.0
cca_rate        = (
    fastq_counters["cca_pass"] / total_al if total_al > 0 else 0.0
)

stats_rows = {
    "sample"              : sample_id,
    "total_pairs"         : fastq_counters["total_input"],
    "both_unmapped"       : fastq_counters["unmapped_written"],
    "both_mapped"         : fastq_counters["mapped_seen"],
    "total_aligned"       : total_al,
    # FIX (read-count reporting): cca_pass/cca_fail below are BAM-level
    # (from Pass 3, one row per alignment) rather than the Pass 2
    # read-name-level counts, since total_aligned/functional_rate are
    # also BAM-level. These can differ from fastq_counters["cca_pass"]/
    # ["cca_fail"] when a read name has multiple alignments.
    "cca_pass"            : bam_counters["cca_pass"],
    "cca_fail"            : bam_counters["cca_fail"],
    # FIX: previously mislabeled — this column held bam_counters["both_pass"]
    # (i.e. it duplicated "both_pass"), not the true anticodon-pass count.
    "anticodon_pass"      : bam_counters["anticodon_pass"],
    "anticodon_fail"      : bam_counters["anticodon_fail"],
    # New: independent breakdown of why a read failed to reach both_pass.
    "fail_cca_only"       : bam_counters["fail_cca_only"],
    "fail_anticodon_only" : bam_counters["fail_anticodon_only"],
    "fail_both"           : bam_counters["fail_both"],
    "both_pass"           : bam_counters["both_pass"],
    "cca_pass_rate"       : round(cca_rate, 4),
    "functional_rate"     : round(functional_rate, 4),
    "unmapped"            : fastq_counters["unmapped_written"],
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
    f"Filter complete. Kept {bam_counters['both_pass']:,} / {total_al:,} aligned reads "
    f"({functional_rate:.1%})."
)
logger.info("CCA+anticodon filter done.")
