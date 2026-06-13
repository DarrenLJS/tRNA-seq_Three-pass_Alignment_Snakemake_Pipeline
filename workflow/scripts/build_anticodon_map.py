import re, logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                    format="%(asctime)s %(message)s")
logger = logging.getLogger()

entries = {}
failed  = []

with open(snakemake.input.fa) as fh:
    for line in fh:
        if not line.startswith(">"):
            continue
        header = line.strip().lstrip(">")
        locus  = header.split()[0]

        parts = locus.split("-")
        if (len(parts) >= 3
                and parts[0].endswith("tRNA")
                and len(parts[2]) == 3
                and re.match(r'^[ACGTU]{3}$', parts[2])):
            entries[locus] = parts[2].replace("U", "T")
            continue

        m = re.search(r'\(([ACGT]{3})\)', header)
        if m:
            entries[locus] = m.group(1)
            continue

        failed.append(locus)
        logger.warning(f"Could not parse anticodon from: {header}")

with open(snakemake.output.tsv, "w") as out:
    out.write("locus\tanticodon\n")
    for locus, ac in sorted(entries.items()):
        out.write(f"{locus}\t{ac}\n")

logger.info(f"Written {len(entries)} locus→anticodon entries to {snakemake.output.tsv}")
if failed:
    logger.warning(f"Failed to parse {len(failed)} entries: {failed[:5]}")
