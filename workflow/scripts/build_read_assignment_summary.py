import pandas as pd, os, re, logging

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
                  "pass1_functional_rate"]:
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

    # ---- QC flags ------------------------------------------------------------
    row["flag_low_cca"] = (
        "WARN" if (row.get("pass1_functional_rate") or 1) < thresh["min_cca_concordance"]
        else "OK"
    )
    row["flag_high_mirna"] = (
        "WARN" if (row.get("pass3_mirna_fraction") or 0) > thresh["max_mirna_fraction"]
        else "OK"
    )

    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(snakemake.output.summary, sep="\t", index=False)
logger.info(f"Written read assignment summary to {snakemake.output.summary}")
logger.info(f"Samples with low CCA rate: {(df['flag_low_cca']=='WARN').sum()}")
logger.info(f"Samples with high miRNA fraction: {(df['flag_high_mirna']=='WARN').sum()}")
