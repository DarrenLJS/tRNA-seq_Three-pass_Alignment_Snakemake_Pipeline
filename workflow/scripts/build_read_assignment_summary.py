import pandas as pd, os, re, glob, logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(message)s")
logger = logging.getLogger()

os.makedirs(os.path.dirname(snakemake.output.summary), exist_ok=True)

rows = []
cl_samples = [
    os.path.basename(p).split(".")[0]           # derive sample IDs from input paths
    for p in snakemake.input.filter_stats
]

scratch   = snakemake.params.scratch
cell_line = snakemake.params.cell_line
thresh    = snakemake.params.qc_thresh
min_pos34_cov = snakemake.params.min_pos34_cov

# ---------------------------------------------------------------------------
# Parse mismatch/ directory once for the whole cell line.
# mim-tRNAseq writes per-position misincorporation data here (rule 03).
# We extract position-34 (wobble) coverage per sample as a QC check:
# sufficient depth at this position is a prerequisite for the Binomial GLM
# wobble-modification inference planned for rule 09.
#
# Expected file: one or more TSVs inside mismatch/ containing at minimum:
#   pos    — canonical tRNA position (1-based; 34 = wobble)
#   cov    — read coverage at this position
#   <id>   — sample identifier column (sample / id / library / name)
#
# NOTE: if the mismatch file format from your mimseq version differs, adjust
# the column name detection in parse_mismatch_coverage() below.
# ---------------------------------------------------------------------------
def parse_mismatch_coverage(mismatch_dir, cl_samples):
    """
    Return dict: sample_id -> {n_isodecoders_pos34_covered, median_pos34_coverage}.
    Returns empty dict on any parse failure; caller fills columns with None.
    """
    candidates = (
        glob.glob(os.path.join(mismatch_dir, "*mismatch*.txt")) +
        glob.glob(os.path.join(mismatch_dir, "*.txt"))
    )
    # Deduplicate while preserving order (mismatch-named files first)
    seen = set()
    candidates = [p for p in candidates if not (p in seen or seen.add(p))]

    if not candidates:
        logger.warning(f"No files found in mismatch dir: {mismatch_dir}")
        return {}

    for fpath in candidates:
        try:
            df = pd.read_csv(fpath, sep="\t")
            df.columns = [c.lower().strip() for c in df.columns]

            if not {"pos", "cov"}.issubset(df.columns):
                logger.info(
                    f"Skipping {os.path.basename(fpath)}: missing pos/cov columns "
                    f"(found: {list(df.columns)})"
                )
                continue

            # Identify the sample-label column (name varies across mimseq versions)
            sample_col = next(
                (c for c in df.columns if c in ("sample", "id", "library", "name")),
                None
            )
            if sample_col is None:
                logger.info(
                    f"Skipping {os.path.basename(fpath)}: no sample identifier column "
                    f"(found: {list(df.columns)})"
                )
                continue

            pos34 = df[df["pos"] == 34].copy()
            if pos34.empty:
                logger.warning(f"No position-34 rows in {os.path.basename(fpath)}")
                continue

            logger.info(
                f"Parsed mismatch data from {os.path.basename(fpath)}: "
                f"{len(pos34)} position-34 rows across "
                f"{pos34[sample_col].nunique()} distinct sample labels"
            )

            result = {}
            for s in cl_samples:
                # mim-tRNAseq labels may retain _val_1 suffix; match if our
                # sample_id is equal to or a prefix of the label in the file.
                mask = pos34[sample_col].apply(
                    lambda x: str(x) == s or str(x).startswith(s)
                )
                s_df = pos34[mask]
                if s_df.empty:
                    logger.warning(f"No position-34 data found for sample {s}")
                    result[s] = {"n_isodecoders_pos34_covered": 0,
                                 "median_pos34_coverage":       0.0}
                else:
                    result[s] = {
                        "n_isodecoders_pos34_covered": int(len(s_df)),
                        "median_pos34_coverage":       round(float(s_df["cov"].median()), 1),
                    }
            return result

        except Exception as e:
            logger.warning(f"Could not parse {os.path.basename(fpath)}: {e}")
            continue

    logger.warning(f"Could not parse any mismatch file in {mismatch_dir}")
    return {}


mismatch_stats = parse_mismatch_coverage(
    snakemake.input.mismatch_dir, cl_samples
)

for s in cl_samples:
    row = {"sample": s, "cell_line": cell_line}

    # ---- Pass 1 filter stats ------------------------------------------------
    fstats_path = f"{scratch}/pass1_filters/{s}/{s}.filter_stats.tsv"
    try:
        fstats = pd.read_csv(fstats_path, sep="\t", index_col=0).squeeze()
        row["pass1_total_aligned"]    = int(fstats.get("total_aligned",   0))
        row["pass1_cca_pass"]         = int(fstats.get("cca_pass",        0))
        row["pass1_anticodon_pass"]   = int(fstats.get("anticodon_pass",  0))
        row["pass1_functional"]       = int(fstats.get("both_pass",       0))
        row["pass1_unmapped"]         = int(fstats.get("unmapped",        0))
        # NEW (feedback item 3): the four independent CCA x anticodon
        # buckets, already computed in cca_anticodon_filter.py's Pass 3
        # (FIX 2026-07-02) but not previously carried forward into this
        # summary -- pass1_cca_pass/pass1_anticodon_pass above are each
        # marginal totals (pass that filter regardless of the other), which
        # cannot be used alone to plot "fail CCA only" vs "fail anticodon
        # only" vs "fail both" vs "pass both". These four columns give the
        # full, mutually-exclusive partition of total_aligned directly.
        row["pass1_both_pass"]           = int(fstats.get("both_pass",           0))
        row["pass1_fail_cca_only"]       = int(fstats.get("fail_cca_only",       0))
        row["pass1_fail_anticodon_only"] = int(fstats.get("fail_anticodon_only", 0))
        row["pass1_fail_both"]           = int(fstats.get("fail_both",           0))
        total = row["pass1_total_aligned"] + row["pass1_unmapped"]
        row["pass1_align_rate"]       = (
            row["pass1_total_aligned"] / total if total > 0 else 0
        )
        row["pass1_functional_rate"]  = (
            row["pass1_functional"] / row["pass1_total_aligned"]
            if row["pass1_total_aligned"] > 0 else 0
        )
    except Exception as e:
        logger.warning(f"Could not parse filter stats for {s}: {e}")
        for k in ["pass1_total_aligned","pass1_cca_pass","pass1_anticodon_pass",
                  "pass1_functional","pass1_unmapped","pass1_align_rate",
                  "pass1_functional_rate","pass1_both_pass","pass1_fail_cca_only",
                  "pass1_fail_anticodon_only","pass1_fail_both"]:
            row[k] = None

    # ---- Pass 2 Bowtie2 stats -----------------------------------------------
    bt2_path = f"{scratch}/pass2_pretRNA/{s}/{s}.bowtie2_stats.txt"
    try:
        with open(bt2_path) as fh:
            content = fh.read()
        m = re.search(r'(\d+) \([\d.]+%\) aligned concordantly exactly 1 time', content)
        row["pass2_mapped_once"] = int(m.group(1)) if m else 0
        m2 = re.search(r'(\d+) \([\d.]+%\) aligned concordantly >1 times', content)
        row["pass2_mapped_multi"] = int(m2.group(1)) if m2 else 0
        row["pass2_mapped_total"] = row["pass2_mapped_once"] + row["pass2_mapped_multi"]
    except Exception as e:
        logger.warning(f"Could not parse Bowtie2 stats for {s}: {e}")
        row["pass2_mapped_total"] = None

    # ---- Pass 3 miRNA fraction -----------------------------------------------
    mir_path = f"{scratch}/pass3_mirna/{s}/{s}_miRNA_fraction.tsv"
    try:
        mir = pd.read_csv(mir_path, sep="\t").iloc[0]
        row["pass3_mirna_reads"]    = int(mir["mirna_reads"])
        row["pass3_mirna_fraction"] = float(mir["mirna_fraction"])
        row["pass3_total_unmapped"] = int(mir["total_pass2_unmapped"])
        row["pass3_unassigned"]     = (
            row["pass3_total_unmapped"] - row["pass3_mirna_reads"]
        )
    except Exception as e:
        logger.warning(f"Could not parse miRNA fraction for {s}: {e}")
        for k in ["pass3_mirna_reads","pass3_mirna_fraction",
                  "pass3_total_unmapped","pass3_unassigned"]:
            row[k] = None

    # ---- Wobble-position 34 mismatch coverage --------------------------------
    mm = mismatch_stats.get(s)
    if mm is not None:
        row["n_isodecoders_pos34_covered"] = mm["n_isodecoders_pos34_covered"]
        row["median_pos34_coverage"]       = mm["median_pos34_coverage"]
    else:
        row["n_isodecoders_pos34_covered"] = None
        row["median_pos34_coverage"]       = None

    # ---- QC flags ------------------------------------------------------------
    row["flag_low_cca"] = (
        "WARN" if (row.get("pass1_functional_rate") or 1) < thresh["min_cca_concordance"]
        else "OK"
    )
    row["flag_high_mirna"] = (
        "WARN" if (row.get("pass3_mirna_fraction") or 0) > thresh["max_mirna_fraction"]
        else "OK"
    )
    med_cov = row.get("median_pos34_coverage")
    row["flag_low_mismatch_cov"] = (
        "WARN" if (med_cov is None or med_cov < min_pos34_cov) else "OK"
    )

    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(snakemake.output.summary, sep="\t", index=False)
logger.info(f"Written read assignment summary to {snakemake.output.summary}")
logger.info(f"Samples with low CCA rate:           {(df['flag_low_cca']=='WARN').sum()}")
logger.info(f"Samples with high miRNA fraction:    {(df['flag_high_mirna']=='WARN').sum()}")
logger.info(f"Samples with low pos-34 mismatch cov:{(df['flag_low_mismatch_cov']=='WARN').sum()}")
