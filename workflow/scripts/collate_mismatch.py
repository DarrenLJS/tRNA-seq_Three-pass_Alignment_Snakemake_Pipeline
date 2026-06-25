"""
workflow/scripts/collate_mismatch.py
=====================================
Snakemake script (called via `script:` directive).
Assembles per-cell-line wobble-position (34) mismatch matrices from
mim-tRNAseq's mismatch/ directory output, for Binomial GLM input.

Directly analogous to collate_counts.py for isodecoder/isoacceptor count
matrices — same input granularity (per cell line), same output structure
(rows = features, cols = samples in manifest order).

The Binomial GLM requires two values per isodecoder per sample:
  - coverage (trials)       → pos34_coverage_matrix.tsv
  - mismatch count (successes) → pos34_mismatch_matrix.tsv

If mim-tRNAseq stores mismatch proportions rather than raw counts,
mismatch count is recovered as: round(proportion * coverage).

Expected mismatch file format (TSV, column names vary across mimseq versions):
  isodecoder/cluster  — feature identifier (row key)
  pos / position      — canonical tRNA position (1-based; filter to 34)
  sample / id         — sample identifier (column key)
  cov / coverage      — read depth at this position (trials)
  mm / mismatch /     — mismatch count OR proportion (successes)
    proportion / rate

Inputs  (snakemake.input):
    mismatch_dir : path to mismatch/ directory produced by rule 03

Outputs (snakemake.output):
    cov_matrix   : pos34_coverage_matrix.tsv
    mm_matrix    : pos34_mismatch_matrix.tsv
"""

import pandas as pd
import os
import glob
import logging

# ---------------------------------------------------------------------------
# Snakemake bindings
# ---------------------------------------------------------------------------
mismatch_dir    = snakemake.input.mismatch_dir
cov_matrix_path = snakemake.output.cov_matrix
mm_matrix_path  = snakemake.output.mm_matrix
cell_line       = snakemake.params.cell_line
manifest_path   = snakemake.params.manifest
log_path        = snakemake.log[0]

os.makedirs(snakemake.params.outdir, exist_ok=True)
os.makedirs(os.path.dirname(log_path), exist_ok=True)

logging.basicConfig(
    filename=log_path, level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger()
logger.info(f"Assembling pos-34 mismatch matrices for cell line: {cell_line}")

# ---------------------------------------------------------------------------
# Load manifest — determine sample order for matrix columns
# ---------------------------------------------------------------------------
manifest   = pd.read_csv(manifest_path, sep="\t", index_col="sample_id")
cl_samples = manifest[manifest["cell_line"] == cell_line].index.tolist()
logger.info(f"Cell line {cell_line}: {len(cl_samples)} samples: {cl_samples}")

# ---------------------------------------------------------------------------
# Locate mismatch data file
# Prefer files with "mismatch" in the name; fall back to any .txt
# ---------------------------------------------------------------------------
candidates = (
    glob.glob(os.path.join(mismatch_dir, "*mismatch*.txt")) +
    glob.glob(os.path.join(mismatch_dir, "*.txt"))
)
seen = set()
candidates = [p for p in candidates if not (p in seen or seen.add(p))]

if not candidates:
    raise FileNotFoundError(f"No files found in mismatch dir: {mismatch_dir}")

logger.info(f"Candidate files: {[os.path.basename(p) for p in candidates]}")

# ---------------------------------------------------------------------------
# Parse mismatch file — flexible column detection
# ---------------------------------------------------------------------------
df         = None
used_file  = None
iso_col    = None
sample_col = None
mm_col     = None

for fpath in candidates:
    try:
        candidate_df = pd.read_csv(fpath, sep="\t")
        candidate_df.columns = [c.lower().strip() for c in candidate_df.columns]

        # Must have position and coverage columns
        if not {"pos", "cov"}.issubset(candidate_df.columns):
            logger.info(
                f"Skipping {os.path.basename(fpath)}: missing pos/cov "
                f"(found: {list(candidate_df.columns)})"
            )
            continue

        # Sample identifier column
        s_col = next(
            (c for c in candidate_df.columns
             if c in ("sample", "id", "library", "name")),
            None
        )
        if s_col is None:
            logger.info(
                f"Skipping {os.path.basename(fpath)}: no sample column "
                f"(found: {list(candidate_df.columns)})"
            )
            continue

        # Isodecoder / feature column
        i_col = next(
            (c for c in candidate_df.columns
             if c in ("isodecoder", "cluster", "trna", "gene", "feature")
             and c != s_col),
            None
        )
        if i_col is None:
            # Fall back: first non-numeric, non-reserved string column
            reserved = {s_col, "pos", "cov"}
            i_col = next(
                (c for c in candidate_df.columns
                 if c not in reserved
                 and candidate_df[c].dtype == object),
                None
            )
        if i_col is None:
            logger.info(
                f"Skipping {os.path.basename(fpath)}: no isodecoder column "
                f"(found: {list(candidate_df.columns)})"
            )
            continue

        # Mismatch count or proportion column
        m_col = next(
            (c for c in candidate_df.columns
             if c in ("mm", "mismatch", "mismatch_count", "count",
                      "proportion", "rate", "mismatch_rate",
                      "mismatch_proportion", "fraction")),
            None
        )
        if m_col is None:
            logger.info(
                f"Skipping {os.path.basename(fpath)}: no mismatch column "
                f"(found: {list(candidate_df.columns)})"
            )
            continue

        df         = candidate_df
        used_file  = fpath
        iso_col    = i_col
        sample_col = s_col
        mm_col     = m_col
        logger.info(
            f"Parsed {os.path.basename(fpath)}: "
            f"iso_col='{iso_col}', sample_col='{sample_col}', "
            f"cov_col='cov', mm_col='{mm_col}'"
        )
        break

    except Exception as e:
        logger.warning(f"Could not parse {os.path.basename(fpath)}: {e}")
        continue

if df is None:
    raise RuntimeError(
        f"Could not parse any mismatch file in {mismatch_dir}. "
        f"Check log for column names found — you may need to update "
        f"the column name lists in collate_mismatch.py."
    )

# ---------------------------------------------------------------------------
# Filter to position 34 (wobble position, 1-based canonical tRNA numbering)
# ---------------------------------------------------------------------------
pos34 = df[df["pos"] == 34].copy()
logger.info(
    f"Position-34 rows: {len(pos34)} across "
    f"{pos34[iso_col].nunique()} isodecoders and "
    f"{pos34[sample_col].nunique()} sample labels"
)

if pos34.empty:
    raise RuntimeError(
        f"No position-34 rows in {os.path.basename(used_file)}. "
        f"Positions present: {sorted(df['pos'].unique()[:20])} "
        f"(showing first 20)"
    )

# ---------------------------------------------------------------------------
# Map mim-tRNAseq sample labels → manifest sample_ids
# Labels may retain _val_1 suffix (e.g. A549_c2_1_val_1 → A549_c2_1)
# ---------------------------------------------------------------------------
def match_sample(label, cl_samples):
    label = str(label)
    for s in cl_samples:
        if label == s or label.startswith(s):
            return s
    return None

pos34 = pos34.copy()
pos34["sample_id"] = pos34[sample_col].apply(lambda x: match_sample(x, cl_samples))

unmatched_mask = pos34["sample_id"].isna()
if unmatched_mask.any():
    logger.warning(
        f"{unmatched_mask.sum()} rows could not be matched to a sample_id — "
        f"labels: {pos34.loc[unmatched_mask, sample_col].unique().tolist()}"
    )
pos34 = pos34[~unmatched_mask].copy()

# ---------------------------------------------------------------------------
# Recover mismatch count if stored as proportion/rate
# Binomial GLM requires integer counts, not rates
# ---------------------------------------------------------------------------
is_proportion = mm_col in (
    "proportion", "rate", "mismatch_rate", "mismatch_proportion", "fraction"
)

if is_proportion:
    logger.info(
        f"'{mm_col}' is a proportion — recovering count as "
        f"round(proportion * coverage)"
    )
    pos34["mm_count"] = (pos34[mm_col] * pos34["cov"]).round().astype(int)
else:
    pos34["mm_count"] = pos34[mm_col].round().astype(int)

# Clamp mismatch count to [0, coverage] — guard against float rounding artefacts
pos34["mm_count"] = pos34["mm_count"].clip(lower=0, upper=pos34["cov"])

# ---------------------------------------------------------------------------
# Pivot to wide matrices: rows = isodecoders, cols = sample_ids
# Use mean to collapse any duplicate (isodecoder, sample) pairs
# ---------------------------------------------------------------------------
cov_wide = (
    pos34.groupby([iso_col, "sample_id"])["cov"]
    .mean()
    .round()
    .astype(int)
    .unstack(fill_value=0)
)
mm_wide = (
    pos34.groupby([iso_col, "sample_id"])["mm_count"]
    .mean()
    .round()
    .astype(int)
    .unstack(fill_value=0)
)

# Ensure all manifest samples are present; fill zeros for any that are absent
for s in cl_samples:
    if s not in cov_wide.columns:
        logger.warning(f"Sample {s} absent from mismatch data — filling with 0")
        cov_wide[s] = 0
        mm_wide[s]  = 0

# Reorder columns to manifest order
cov_wide = cov_wide[cl_samples]
mm_wide  = mm_wide[cl_samples]

logger.info(
    f"Coverage matrix : {cov_wide.shape[0]} isodecoders × "
    f"{cov_wide.shape[1]} samples"
)
logger.info(
    f"Mismatch matrix : {mm_wide.shape[0]} isodecoders × "
    f"{mm_wide.shape[1]} samples"
)

# ---------------------------------------------------------------------------
# Validation: row indices must match between the two matrices
# ---------------------------------------------------------------------------
assert list(cov_wide.index) == list(mm_wide.index), (
    "Row mismatch between coverage and mismatch matrices — "
    "check for duplicate (isodecoder, sample) pairs in the mismatch file"
)
assert list(cov_wide.columns) == cl_samples, \
    "Column order mismatch between coverage matrix and manifest"
assert list(mm_wide.columns) == cl_samples, \
    "Column order mismatch between mismatch matrix and manifest"

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
cov_wide.to_csv(cov_matrix_path, sep="\t")
mm_wide.to_csv(mm_matrix_path,   sep="\t")

logger.info(f"Written coverage matrix → {cov_matrix_path}")
logger.info(f"Written mismatch matrix → {mm_matrix_path}")
logger.info("collate_mismatch.py complete.")
