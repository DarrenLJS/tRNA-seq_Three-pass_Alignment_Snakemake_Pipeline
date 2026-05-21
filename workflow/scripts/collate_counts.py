"""
workflow/scripts/collate_counts.py
====================================
Snakemake script (called via `script:` directive).
Builds DESeq2-ready count matrices and coldata from mim-tRNAseq output.

mim-tRNAseq outputs a single Isodecoder_counts.txt and Isoacceptor_counts.txt
per cell line, with columns = samples and rows = tRNA features.
This script:
  1. Reads both count files
  2. Renames columns to match sample_ids from the manifest
  3. Writes the matrices in the format expected by DESeq2 (rows=features, cols=samples)
  4. Writes coldata.tsv: sample metadata (condition, timepoint, replicate)
     matching the column order of the count matrix

DESeq2 design formula (for later use in R):
    ~ condition + timepoint + condition:timepoint
or for paired time-course LRT:
    full: ~ condition * timepoint
    reduced: ~ timepoint

Inputs  (snakemake.input):
    iso_counts   : mim-tRNAseq Isodecoder_counts.txt
    isoa_counts  : mim-tRNAseq Isoacceptor_counts.txt

Outputs (snakemake.output):
    iso_matrix   : isodecoder count matrix TSV
    isoa_matrix  : isoacceptor count matrix TSV
    coldata      : sample metadata TSV
"""

import pandas as pd
import os
import logging

# ---------------------------------------------------------------------------
# Snakemake bindings
# ---------------------------------------------------------------------------
iso_counts_path   = snakemake.input.iso_counts
isoa_counts_path  = snakemake.input.isoa_counts
iso_matrix_path   = snakemake.output.iso_matrix
isoa_matrix_path  = snakemake.output.isoa_matrix
coldata_path      = snakemake.output.coldata
cell_line         = snakemake.params.cell_line
log_path          = snakemake.log[0]

# The global manifest is available in the Snakemake rule's scope
# We re-read it here since scripts run in a separate Python process
manifest_path = snakemake.config["manifest"]

os.makedirs(snakemake.params.outdir, exist_ok=True)
os.makedirs(os.path.dirname(log_path), exist_ok=True)

logging.basicConfig(
    filename=log_path, level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s"
)
logger = logging.getLogger()
logger.info(f"Collating counts for cell line: {cell_line}")

# ---------------------------------------------------------------------------
# Load manifest for this cell line
# ---------------------------------------------------------------------------
manifest = pd.read_csv(manifest_path, sep="\t", index_col="sample_id")
cl_manifest = manifest[manifest["cell_line"] == cell_line].copy()
cl_samples  = cl_manifest.index.tolist()
logger.info(f"Cell line {cell_line}: {len(cl_samples)} samples")
logger.info(f"Samples: {cl_samples}")


def load_and_rename(counts_path, cl_samples):
    """
    Load a mim-tRNAseq counts file, rename columns to match sample_ids,
    and return the DataFrame.

    mim-tRNAseq column names are based on input FASTQ basenames
    (stripping _val_1.fq.gz). Sample columns should match {sample}_val_1
    or similar. We match by sample_id (flexible to handle naming variants).
    """
    df = pd.read_csv(counts_path, sep="\t", index_col=0)
    logger.info(f"Loaded {counts_path}: {df.shape[0]} features × {df.shape[1]} samples")
    logger.info(f"Column names in file: {df.columns.tolist()}")

    # Build column rename map: file_column → sample_id
    # mim-tRNAseq strips extensions from input filenames to form column headers.
    # Expected pattern: {sample_id}_val_1 or {sample_id}
    rename_map = {}
    for col in df.columns:
        # Try direct match first
        if col in cl_samples:
            rename_map[col] = col
            continue
        # Try stripping _val_1 suffix
        stripped = col.replace("_val_1", "").replace("_val_2", "")
        if stripped in cl_samples:
            rename_map[col] = stripped
            continue
        # Warn if no match
        logger.warning(f"Could not match column '{col}' to any sample_id")

    df = df.rename(columns=rename_map)

    # Reorder columns to match manifest order
    matched_cols = [s for s in cl_samples if s in df.columns]
    missing = [s for s in cl_samples if s not in df.columns]
    if missing:
        logger.error(f"Missing sample columns after rename: {missing}")

    df = df[matched_cols]

    # Drop rows where all counts are 0
    all_zero = (df == 0).all(axis=1)
    n_zero   = all_zero.sum()
    if n_zero > 0:
        logger.info(f"Dropping {n_zero} all-zero features")
    df = df[~all_zero]

    logger.info(f"Final matrix: {df.shape[0]} features × {df.shape[1]} samples")
    return df


# ---------------------------------------------------------------------------
# Load and process count tables
# ---------------------------------------------------------------------------
iso_df  = load_and_rename(iso_counts_path,  cl_samples)
isoa_df = load_and_rename(isoa_counts_path, cl_samples)

# Save count matrices
iso_df.to_csv(iso_matrix_path,  sep="\t")
isoa_df.to_csv(isoa_matrix_path, sep="\t")
logger.info(f"Written isodecoder matrix  → {iso_matrix_path}")
logger.info(f"Written isoacceptor matrix → {isoa_matrix_path}")

# ---------------------------------------------------------------------------
# Build coldata
# ---------------------------------------------------------------------------
# Columns in coldata must match columns in the count matrix (same order).
coldata = cl_manifest.loc[iso_df.columns, ["condition", "timepoint", "replicate"]].copy()
coldata.index.name = "sample_id"

# DESeq2 requires factors; encode condition as factor-friendly strings
coldata["condition"] = coldata["condition"].astype(str)
coldata["timepoint"] = coldata["timepoint"].astype(str)
coldata["replicate"] = coldata["replicate"].astype(str)

# Add a combined group label for simple pairwise contrasts
coldata["group"] = coldata["condition"] + "_" + coldata["timepoint"] + "h"

coldata.to_csv(coldata_path, sep="\t")
logger.info(f"Written coldata → {coldata_path}")
logger.info(f"\n{coldata.to_string()}")

# ---------------------------------------------------------------------------
# Validation: check column order consistency between matrix and coldata
# ---------------------------------------------------------------------------
assert list(iso_df.columns)  == list(coldata.index), \
    "Column order mismatch between isodecoder matrix and coldata!"
assert list(isoa_df.columns) == list(coldata.index), \
    "Column order mismatch between isoacceptor matrix and coldata!"

logger.info("Column order validation passed.")
logger.info("collate_counts.py complete.")
