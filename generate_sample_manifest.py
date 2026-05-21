#!/usr/bin/python3

# generate_sample_manifest.py
import os, re, csv

CELL_LINES = {
    "A549": "/scratch/s2906787/A549/raw_files",
    "THP1": "/scratch/s2906787/THP1/raw_files",
}

# A549 pattern: lane and read encoded in filename
# THP1 pattern: one lane, _L1_1.fq.gz / _L1_2.fq.gz

rows = []

for cell_line, base in CELL_LINES.items():
    for sample_dir in sorted(os.listdir(base)):
        path = os.path.join(base, sample_dir)
        if not os.path.isdir(path):
            continue

        # Parse sample name
        name = sample_dir.replace("sg5_", "")   # normalise THP1 prefix
        m = re.match(r"([cp])(\d+)_(\d+)", name)
        if not m:
            continue
        condition = "control" if m.group(1) == "c" else "polyIC"
        timepoint = int(m.group(2))
        replicate = int(m.group(3))

        files = sorted([
            os.path.join(path, f)
            for f in os.listdir(path)
            if f.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz"))
        ])

        # Group by R1/R2
        r1 = [f for f in files if "_1.fastq" in f or "_1.fq" in f
                                or "_1.fastq.gz" in f or "_1.fq.gz" in f]
        r2 = [f for f in files if "_2.fastq" in f or "_2.fq" in f
                                or "_2.fastq.gz" in f or "_2.fq.gz" in f]

        rows.append({
            "sample_id":  f"{cell_line}_{name}",
            "cell_line":  cell_line,
            "condition":  condition,
            "timepoint":  timepoint,
            "replicate":  replicate,
            "R1_files":   ";".join(r1),   # semicolon-separated for multi-lane
            "R2_files":   ";".join(r2),
            "n_R1_files": len(r1),
        })

with open("sample_manifest.tsv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Written {len(rows)} samples to sample_manifest.tsv")