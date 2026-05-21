# =============================================================================
# workflow/rules/00_reference_prep.smk
# Build all derived reference files needed by downstream rules.
#
# PREREQUISITES (download manually before first run):
#   1. GRCh38 primary assembly FASTA
#      wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/\
#           Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#      gunzip && samtools faidx
#
#   2. GtRNAdb Hsapi38 files (BED + FASTA + introns)
#      https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/
#      Download: hg38-tRNAs.bed, hg38-tRNAs-detailed.out,
#                hg38-tRNA-introns.bed, hg38-tRNAs.fa
#
#   3. miRBase release 22.1
#      wget https://www.mirbase.org/download/mature.fa
#      wget https://www.mirbase.org/download/hairpin.fa
# =============================================================================

REF = config["references"]


# ---------------------------------------------------------------------------
# 00a: Build flanked pre-tRNA FASTA (intron-retaining form)
#      Extends each tRNA locus by 50 nt on each side using bedtools slop,
#      then extracts sequence with bedtools getfasta.
# ---------------------------------------------------------------------------
rule build_pretRNA_fasta:
    """
    Construct the pre-tRNA reference for Pass 2 Bowtie2 alignment.
    Each entry is a tRNA locus extended ±50 nt to capture intronic and
    flanking reads (as per proposal Section 3.3, Pass 2).
    """
    input:
        bed    = REF["gtrndb_bed"],
        fasta  = REF["genome_fasta"],
        fai    = REF["genome_fai"],
    output:
        flanked_bed  = REF["gtrndb_bed"].replace(".bed", "_flanked50.bed"),
        pretRNA_fa   = REF["pretRNA_fasta"],
    log:
        f"{SCRATCH}/logs/00_build_pretRNA_fasta.log",
    benchmark:
        f"{SCRATCH}/benchmarks/00_build_pretRNA_fasta.tsv",
    shell:
        r"""
        set -euo pipefail
        exec &> {log}

        echo "[$(date)] Slopping tRNA BED ±50 nt..."
        bedtools slop \
            -i   {input.bed} \
            -g   {input.fai} \
            -b   50 \
            > {output.flanked_bed}

        echo "[$(date)] Extracting pre-tRNA FASTA (intron-retaining)..."
        bedtools getfasta \
            -fi   {input.fasta} \
            -bed  {output.flanked_bed} \
            -name \
            -s \
            > {output.pretRNA_fa}

        echo "[$(date)] Done. $(grep -c '^>' {output.pretRNA_fa}) sequences written."
        """


# ---------------------------------------------------------------------------
# 00b: Build spliced pre-tRNA FASTA (for intron-containing loci)
#      Removes intronic intervals then extracts exon sequence.
#      Both forms (00a and 00b) are cat'd into the final Bowtie2 reference.
# ---------------------------------------------------------------------------
rule build_pretRNA_spliced_fasta:
    """
    Build a spliced pre-tRNA FASTA by subtracting intron intervals from
    the flanked loci, producing exon-only sequences for intron-containing genes.
    Merged with the intron-retaining form so Bowtie2 can capture both.
    """
    input:
        flanked_bed  = REF["gtrndb_bed"].replace(".bed", "_flanked50.bed"),
        introns      = REF["gtrndb_introns"],
        fasta        = REF["genome_fasta"],
    output:
        spliced_bed  = REF["gtrndb_bed"].replace(".bed", "_flanked50_spliced.bed"),
        spliced_fa   = REF["pretRNA_fasta_spliced"],
    log:
        f"{SCRATCH}/logs/00_build_pretRNA_spliced_fasta.log",
    benchmark:
        f"{SCRATCH}/benchmarks/00_build_pretRNA_spliced_fasta.tsv",
    shell:
        r"""
        set -euo pipefail
        exec &> {log}

        echo "[$(date)] Subtracting introns from flanked pre-tRNA loci..."
        bedtools subtract \
            -a {input.flanked_bed} \
            -b {input.introns} \
            > {output.spliced_bed}

        echo "[$(date)] Extracting spliced pre-tRNA FASTA..."
        bedtools getfasta \
            -fi   {input.fasta} \
            -bed  {output.spliced_bed} \
            -name \
            -s \
            > {output.spliced_fa}

        echo "[$(date)] Done. $(grep -c '^>' {output.spliced_fa}) sequences written."
        """


# ---------------------------------------------------------------------------
# 00c: Build combined (intron-retaining + spliced) Bowtie2 index
# ---------------------------------------------------------------------------
rule build_bowtie2_pretRNA_index:
    """
    Concatenate intron-retaining and spliced pre-tRNA FASTAs, then build
    a Bowtie2 index from the combined reference.
    """
    input:
        fa_intron  = REF["pretRNA_fasta"],
        fa_spliced = REF["pretRNA_fasta_spliced"],
    output:
        combined   = REF["pretRNA_fasta"].replace(".fa", "_combined.fa"),
        index_done = REF["pretRNA_index"] + ".1.bt2",
    params:
        index_prefix = REF["pretRNA_index"],
        threads      = config["bowtie2"]["threads"],
    log:
        f"{SCRATCH}/logs/00_build_bowtie2_pretRNA_index.log",
    benchmark:
        f"{SCRATCH}/benchmarks/00_build_bowtie2_pretRNA_index.tsv",
    shell:
        r"""
        set -euo pipefail
        exec &> {log}

        echo "[$(date)] Combining intron-retaining and spliced pre-tRNA FASTAs..."
        cat {input.fa_intron} {input.fa_spliced} > {output.combined}

        echo "[$(date)] Building Bowtie2 index at {params.index_prefix}..."
        mkdir -p $(dirname {params.index_prefix})
        bowtie2-build \
            --threads {params.threads} \
            {output.combined} \
            {params.index_prefix}

        echo "[$(date)] Bowtie2 index built."
        """


# ---------------------------------------------------------------------------
# 00d: Build anticodon lookup table from GtRNAdb FASTA headers
#      Output TSV: locus_name <TAB> anticodon (3-char, e.g. AGC)
#      Used by the pysam anticodon concordance filter (Rule 04).
# ---------------------------------------------------------------------------
rule build_anticodon_map:
    """
    Parse GtRNAdb FASTA headers to extract the locus → anticodon mapping.

    GtRNAdb 2.0 NEW name format (what we have):
      >tRNA-Tyr-GTA-1-1  (tRNA-{AminoAcid}-{Anticodon}-{Family}-{Copy})
      The anticodon is the third dash-separated field.

    Also handles old format as fallback:
      >chr1.tRNA1-AlaAGC chr1:... AlaAGC (AGC)

    Output TSV: locus_name <TAB> anticodon
    e.g.: tRNA-Tyr-GTA-1-1    GTA
    """
    input:
        fa = REF["gtrndb_fasta"],
    output:
        tsv = REF["anticodon_map"],
    log:
        f"{SCRATCH}/logs/00_build_anticodon_map.log",
    run:
        import re, logging
        logging.basicConfig(filename=log[0], level=logging.INFO,
                            format="%(asctime)s %(message)s")
        logger = logging.getLogger()

        entries = {}
        failed  = []

        with open(input.fa) as fh:
            for line in fh:
                if not line.startswith(">"):
                    continue
                header = line.strip().lstrip(">")
                locus  = header.split()[0]   # first whitespace-delimited field

                # ── New GtRNAdb 2.0 format: tRNA-Tyr-GTA-1-1 ──────────────
                # Split on dash: ['tRNA', 'Tyr', 'GTA', '1', '1']
                # Anticodon is index 2, always 3 characters
                parts = locus.split("-")
                if (len(parts) >= 3
                        and parts[0] == "tRNA"
                        and len(parts[2]) == 3
                        and re.match(r'^[ACGTU]{3}$', parts[2])):
                    entries[locus] = parts[2].replace("U", "T")
                    continue

                # ── Old GtRNAdb format fallback: (AGC) at end of header ────
                m = re.search(r'\(([ACGT]{3})\)\s*$', header)
                if m:
                    entries[locus] = m.group(1)
                    continue

                # ── Could not parse ────────────────────────────────────────
                failed.append(locus)
                logger.warning(f"Could not parse anticodon from: {header}")

        with open(output.tsv, "w") as out:
            out.write("locus\tanticodon\n")
            for locus, ac in sorted(entries.items()):
                out.write(f"{locus}\t{ac}\n")

        logger.info(f"Written {len(entries)} locus→anticodon entries to {output.tsv}")
        if failed:
            logger.warning(f"Failed to parse {len(failed)} entries: {failed[:5]}")


# ---------------------------------------------------------------------------
# 00e: Filter miRBase FASTA to human (hsa) entries only
# ---------------------------------------------------------------------------
rule filter_mirbase_human:
    """
    Subset mature.fa and hairpin.fa to human (hsa) sequences only.
    miRDeep2 requires a genome-matched reference; using only hsa avoids
    spurious cross-species alignments.
    """
    input:
        mature  = REF["mirbase_mature_all"],
        hairpin = REF["mirbase_hairpin_all"],
    output:
        mature_hsa  = REF["mirbase_mature_hsa"],
        hairpin_hsa = REF["mirbase_hairpin_hsa"],
    log:
        f"{SCRATCH}/logs/00_filter_mirbase_human.log",
    run:
        import logging
        logging.basicConfig(filename=log[0], level=logging.INFO,
                            format="%(asctime)s %(message)s")

        def filter_fasta(src, dst, prefix="hsa-"):
            kept = 0
            with open(src) as fin, open(dst, "w") as fout:
                write = False
                for line in fin:
                    if line.startswith(">"):
                        write = line[1:].startswith(prefix)
                        if write:
                            kept += 1
                    if write:
                        fout.write(line)
            return kept

        n_mat = filter_fasta(input.mature,  output.mature_hsa)
        n_hp  = filter_fasta(input.hairpin, output.hairpin_hsa)
        logging.info(f"Kept {n_mat} mature hsa miRNAs, {n_hp} hairpin hsa miRNAs.")


# ---------------------------------------------------------------------------
# 00f: Build TRAX database from GtRNAdb files
#      TRAX (tRAX) requires a pre-built database for tRF quantification.
# ---------------------------------------------------------------------------
rule build_trax_db:
    """
    Construct the TRAX database from GtRNAdb files.
    Requires: trax (tRAX) installed in the trax conda env.

    TODO: Verify the exact trax database-build command for your tRAX version.
          Typical usage: maketraxdb.py --trnaout hg38-tRNAs-detailed.out
                                       --genome  GRCh38.fa
                                       --name    hsapi38
    """
    input:
        tRNA_out = REF["gtrndb_detailed_out"],
        genome   = REF["genome_fasta"],
    output:
        db_flag  = directory(REF["trax_ref_dir"]),
    params:
        db_name  = "hsapi38",
        outdir   = REF["trax_ref_dir"],
    log:
        f"{SCRATCH}/logs/00_build_trax_db.log",
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        mkdir -p {params.outdir}

        echo "[$(date)] Building TRAX database..."
        maketraxdb.py \
            --trnaout {input.tRNA_out} \
            --genome  {input.genome} \
            --name    {params.db_name} \
            --out     {params.outdir}

        echo "[$(date)] TRAX database built at {params.outdir}."
        """
