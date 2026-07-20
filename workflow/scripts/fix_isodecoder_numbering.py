"""
fix_isodecoder_numbering.py

GtRNAdb's isodecoder numbering does not guarantee that loci sharing 100%
identical mature sequence share the same isodecoder number. Two loci with
different isodecoder numbers (e.g. Val-CAC-1-* and Val-CAC-5-1) can turn
out to have byte-identical mature sequence.

mim-tRNAseq clusters reads with usearch at --cluster-id 0.97, which
clusters purely on sequence identity and has no knowledge of the isodecoder
number encoded in the locus name. If usearch merges all members of one
isodecoder entirely into another isodecoder's cluster, mim-tRNAseq's
isodecoder-name-based bookkeeping (writeIsodecoderInfo in splitClusters.py)
looks up the now-empty isodecoder and crashes:

    ValueError: min() arg is an empty sequence

(and a related IndexError further downstream in mmQuant.py's multiprocessing
worker, from the same class of name/cluster-membership mismatch).

This script renumbers isodecoders by exact sequence identity within each
(isotype, anticodon) isoacceptor family, so that mim-tRNAseq's name-based
bookkeeping can never disagree with usearch's sequence-based clustering
again. Only the isodecoder and copy number fields in the FASTA header are
changed; the tRNAscan-SE ID, genomic coordinates, strand, and score are
left untouched. Locus order in the file is preserved.

Consumed as a Snakemake `script:` directive — see rule
fix_isodecoder_numbering in 00_reference_prep.smk.
"""

import re
import logging
from collections import defaultdict, OrderedDict

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
                     format="%(asctime)s %(message)s")
logger = logging.getLogger()

# Matches: <prefix>_tRNA-<isotype>-<anticodon>-<isodecoder>-<copy>
# e.g. "Homo_sapiens_tRNA-Ala-AGC-1-1" ->
#      prefix="Homo_sapiens_tRNA", iso="Ala", anti="AGC", isodec=1, copy=1
NAME_RE = re.compile(
    r'^(?P<prefix>.+_tRNA)-(?P<iso>[^-]+)-(?P<anti>[^-]+)-(?P<isodec>\d+)-(?P<copy>\d+)$'
)


def parse_fasta(path):
    """Parse a FASTA file into (name, header_rest, sequence) triples,
    preserving file order. header_rest is everything on the header line
    after the first whitespace-delimited token (name)."""
    records = []
    name, header_rest, seq_lines = None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, header_rest, "".join(seq_lines)))
                fields = line[1:].split(None, 1)
                name = fields[0]
                header_rest = fields[1] if len(fields) > 1 else ""
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if name is not None:
            records.append((name, header_rest, "".join(seq_lines)))
    return records


def wrap_seq(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def renumber(records):
    """
    Returns:
        output_by_idx   : dict idx -> (new_name, header_rest, seq), covers
                           every input record (unparsed records pass through
                           unchanged).
        mapping_rows     : list of (old_name, new_name, changed_bool)
        collision_rows   : list of (family, old_name, new_name, group_size,
                            old_isodecoders_merged) for every record whose
                            sequence group spanned more than one original
                            isodecoder number.
        unparsed         : list of (name, header_rest, seq) that didn't match
                            the expected 5-field naming convention.
    """
    parsed = []
    unparsed = []
    name_to_idx = {}
    for idx, (name, header_rest, seq) in enumerate(records):
        name_to_idx[name] = idx
        m = NAME_RE.match(name)
        if not m:
            unparsed.append((name, header_rest, seq))
            continue
        parsed.append({
            "idx": idx, "name": name, "header_rest": header_rest, "seq": seq,
            "prefix": m.group("prefix"), "iso": m.group("iso"), "anti": m.group("anti"),
            "old_isodec": int(m.group("isodec")), "old_copy": int(m.group("copy")),
        })

    families = defaultdict(list)
    for rec in parsed:
        families[(rec["prefix"], rec["iso"], rec["anti"])].append(rec)

    mapping_rows = []
    collision_rows = []
    output_by_idx = {}

    for fam_key, fam_recs in families.items():
        # Group by exact (case-insensitive) sequence identity within the family.
        seq_groups = OrderedDict()
        for rec in fam_recs:
            seq_groups.setdefault(rec["seq"].upper(), []).append(rec)

        def group_sort_key(item):
            _, members = item
            size = len(members)
            min_old_isodec = min(r["old_isodec"] for r in members)
            min_idx = min(r["idx"] for r in members)
            # Largest group first (mirrors GtRNAdb convention of isodecoder 1
            # typically being the most common class); ties broken by lowest
            # original isodecoder number, then by first file appearance —
            # both fully deterministic.
            return (-size, min_old_isodec, min_idx)

        ordered_groups = sorted(seq_groups.items(), key=group_sort_key)

        for new_isodec, (seq, members) in enumerate(ordered_groups, start=1):
            members_sorted = sorted(members, key=lambda r: r["idx"])
            old_isodecs_in_group = sorted(set(r["old_isodec"] for r in members_sorted))
            merged = len(old_isodecs_in_group) > 1

            for new_copy, rec in enumerate(members_sorted, start=1):
                new_name = f'{rec["prefix"]}-{rec["iso"]}-{rec["anti"]}-{new_isodec}-{new_copy}'
                output_by_idx[rec["idx"]] = (new_name, rec["header_rest"], rec["seq"])
                changed = new_name != rec["name"]
                mapping_rows.append((rec["name"], new_name, changed))
                if merged:
                    collision_rows.append((
                        f'{rec["prefix"]}-{rec["iso"]}-{rec["anti"]}',
                        rec["name"], new_name, len(members_sorted),
                        ",".join(str(x) for x in old_isodecs_in_group),
                    ))

    for name, header_rest, seq in unparsed:
        idx = name_to_idx[name]
        output_by_idx[idx] = (name, header_rest, seq)
        mapping_rows.append((name, name, False))

    return output_by_idx, mapping_rows, collision_rows, unparsed


def main():
    in_fa = snakemake.input.raw_fa
    out_fa = snakemake.output.fixed_fa
    out_mapping = snakemake.output.mapping
    out_report = snakemake.output.report

    logger.info(f"Reading raw reference: {in_fa}")
    records = parse_fasta(in_fa)
    logger.info(f"Parsed {len(records)} FASTA records.")

    output_by_idx, mapping_rows, collision_rows, unparsed = renumber(records)

    if unparsed:
        logger.warning(
            f"{len(unparsed)} record(s) did not match the expected "
            f"<prefix>_tRNA-<iso>-<anticodon>-<isodecoder>-<copy> naming "
            f"convention and were passed through unchanged: "
            f"{[n for n, _, _ in unparsed]}"
        )

    logger.info(f"Writing renumbered reference: {out_fa}")
    with open(out_fa, "w") as out:
        for idx in sorted(output_by_idx):
            new_name, header_rest, seq = output_by_idx[idx]
            header = f">{new_name} {header_rest}" if header_rest else f">{new_name}"
            out.write(header + "\n")
            out.write(wrap_seq(seq) + "\n")

    n_changed = sum(1 for _, _, changed in mapping_rows if changed)
    logger.info(f"Writing name mapping ({len(mapping_rows)} records, "
                f"{n_changed} renamed): {out_mapping}")
    with open(out_mapping, "w") as out:
        out.write("old_name\tnew_name\tchanged\n")
        for old, new, changed in sorted(mapping_rows):
            out.write(f"{old}\t{new}\t{changed}\n")

    n_families_merged = len(set(row[0] for row in collision_rows))
    logger.info(f"Writing collision report ({len(collision_rows)} records "
                f"across {n_families_merged} isoacceptor famil"
                f"{'y' if n_families_merged == 1 else 'ies'} with a cross-"
                f"isodecoder sequence collision): {out_report}")
    with open(out_report, "w") as out:
        out.write("family\told_name\tnew_name\tgroup_size\told_isodecoders_merged\n")
        for row in sorted(collision_rows):
            out.write("\t".join(str(x) for x in row) + "\n")

    if collision_rows:
        families = sorted(set(row[0] for row in collision_rows))
        logger.warning(
            f"Cross-isodecoder sequence collisions found and corrected in "
            f"{len(families)} isoacceptor famil"
            f"{'y' if len(families) == 1 else 'ies'}: {families}. "
            f"See {out_report} for full detail."
        )
    else:
        logger.info("No cross-isodecoder sequence collisions found.")

    logger.info("Done.")


main()
