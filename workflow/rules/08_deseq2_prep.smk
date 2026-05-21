# =============================================================================
# workflow/rules/08_deseq2_prep.smk
# Collate mim-tRNAseq count tables into per-cell-line count matrices and
# write the corresponding coldata (sample metadata) TSVs for DESeq2/edgeR.
#
# Also assembles a per-cell-line read assignment summary (QC metric 2 from
# the proposal: proportional read assignment across all categories).
#
# Outputs per cell line:
#   isodecoder_counts_matrix.tsv   — rows=isodecoders, cols=samples
#   isoacceptor_counts_matrix.tsv  — rows=isoacceptors, cols=samples
#   coldata.tsv                    — sample metadata for DESeq2
#   read_assignment_summary.tsv    — read proportions per sample (QC)
# =============================================================================

rule build_deseq2_inputs:
    """
    Collate per-cell-line count matrices and coldata for DESeq2.
    Calls workflow/scripts/collate_counts.py.
    """
    input:
        iso_counts   = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isodecoder_counts.txt",
        isoa_counts  = f"{SCRATCH}/pass1_mimtrnaseq/{{cell_line}}/counts/Isoacceptor_counts.txt",
    output:
        iso_matrix   = f"{SCRATCH}/deseq2_input/{{cell_line}}/isodecoder_counts_matrix.tsv",
        isoa_matrix  = f"{SCRATCH}/deseq2_input/{{cell_line}}/isoacceptor_counts_matrix.tsv",
        coldata      = f"{SCRATCH}/deseq2_input/{{cell_line}}/coldata.tsv",
    params:
        outdir       = f"{SCRATCH}/deseq2_input/{{cell_line}}",
        cell_line    = "{cell_line}",
    log:
        f"{SCRATCH}/logs/08_deseq2_prep/{{cell_line}}.log",
    script:
        "../scripts/collate_counts.py"


rule build_read_assignment_summary:
    """
    Compile a per-sample read-assignment table across all 3 passes and
    miRNA screening. This is QC metric 2 from the proposal: proportional
    read assignment across mature tRNA / pre-tRNA / miRNA / unassigned.

    Inputs:
      - Trim Galore reports        → total trimmed read counts
      - mim-tRNAseq align_done     → pass 1 mapped counts (from BAM flagstat)
      - Bowtie2 stats files        → pass 2 mapped counts
      - miRDeep2 qc_summary files  → pass 3 miRNA counts
    """
    input:
        filter_stats  = lambda wildcards: [
            f"{SCRATCH}/pass1_filters/{s}/{s}.filter_stats.tsv"
            for s in samples_for(wildcards.cell_line)
        ],
        bowtie2_stats = lambda wildcards: [
            f"{SCRATCH}/pass2_pretRNA/{s}/{s}.bowtie2_stats.txt"
            for s in samples_for(wildcards.cell_line)
        ],
        mirna_frac    = lambda wildcards: [
            f"{SCRATCH}/pass3_mirna/{s}/{s}_miRNA_fraction.tsv"
            for s in samples_for(wildcards.cell_line)
        ],
    output:
        summary = f"{SCRATCH}/qc/read_assignment/{{cell_line}}_read_assignment_summary.tsv",
    params:
        cell_line = "{cell_line}",
        scratch   = SCRATCH,
        qc_thresh = config["qc_thresholds"],
    log:
        f"{SCRATCH}/logs/08_read_assignment/{{cell_line}}.log",
    run:
        import pandas as pd, os, re, logging

        logging.basicConfig(filename=log[0], level=logging.INFO,
                            format="%(asctime)s %(message)s")
        logger = logging.getLogger()

        os.makedirs(os.path.dirname(output.summary), exist_ok=True)

        rows = []
        cl_samples = samples_for(params.cell_line)

        for s in cl_samples:
            row = {"sample": s, "cell_line": params.cell_line}

            # ---- Pass 1 filter stats ----------------------------------------
            fstats_path = f"{params.scratch}/pass1_filters/{s}/{s}.filter_stats.tsv"
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

            # ---- Pass 2 Bowtie2 stats ----------------------------------------
            bt2_path = f"{params.scratch}/pass2_pretRNA/{s}/{s}.bowtie2_stats.txt"
            try:
                with open(bt2_path) as fh:
                    content = fh.read()
                # Parse concordantly aligned once or multiple times
                m = re.search(
                    r'(\d+) \([\d.]+%\) aligned concordantly exactly 1 time', content
                )
                row["pass2_mapped_once"] = int(m.group(1)) if m else 0
                m2 = re.search(
                    r'(\d+) \([\d.]+%\) aligned concordantly >1 times', content
                )
                row["pass2_mapped_multi"] = int(m2.group(1)) if m2 else 0
                row["pass2_mapped_total"] = row["pass2_mapped_once"] + row["pass2_mapped_multi"]
            except Exception as e:
                logger.warning(f"Could not parse Bowtie2 stats for {s}: {e}")
                row["pass2_mapped_total"] = None

            # ---- Pass 3 miRNA fraction ----------------------------------------
            mir_path = f"{params.scratch}/pass3_mirna/{s}/{s}_miRNA_fraction.tsv"
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

            # ---- QC flags ────────────────────────────────────────────────────
            thresh = params.qc_thresh
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
        df.to_csv(output.summary, sep="\t", index=False)
        logger.info(f"Written read assignment summary to {output.summary}")
        logger.info(f"Samples with low CCA rate: "
                    f"{(df['flag_low_cca']=='WARN').sum()}")
        logger.info(f"Samples with high miRNA fraction: "
                    f"{(df['flag_high_mirna']=='WARN').sum()}")
