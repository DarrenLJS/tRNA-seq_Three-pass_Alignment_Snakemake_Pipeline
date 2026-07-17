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
#
# Conda environments used
# -----------------------
#   Rules 00a–00e  →  ../../envs/environment.yaml  (bedtools, bowtie2, samtools)
#   Rule  00f      →  ../../envs/trax_env.yaml      (maketraxdb.py / tRAX)
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
    resources:
        runtime = 60,
    conda:
        "../../envs/environment.yaml"
    shell:
        r"""
        set -euo pipefail
        exec &> {log}

        # NOTE: GRCh38.primary_assembly.genome.fa is from Ensembl and uses
        # chromosome names WITHOUT a chr prefix (e.g. "1", "2").
        # hg38-tRNAs_nochr.bed from GtRNAdb also lacks the chr prefix, so the
        # primary chromosomes already match.  However, GtRNAdb scaffold names
        # (e.g. "1_KI270713v1_random") differ from Ensembl's naming
        # (e.g. "KI270713.1"), causing bedtools slop to error.
        # We pre-filter the BED to only rows whose chromosome appears in the
        # FAI.  These scaffold tRNAs are a tiny minority and dropping them is
        # standard practice for this reference combination.
        REFDIR=$(dirname {output.flanked_bed})
        mkdir -p "$REFDIR"

        echo "[$(date)] Filtering BED to chromosomes present in FAI..."
        awk 'NR==FNR {{chroms[$1]=1; next}} ($1 in chroms)' \
            {input.fai} {input.bed} \
            > "$REFDIR/tRNAs_fai_filtered.bed"
        N_IN=$(wc -l < {input.bed})
        N_OUT=$(wc -l < "$REFDIR/tRNAs_fai_filtered.bed")
        echo "[$(date)] Kept $N_OUT of $N_IN BED entries after chromosome filter."

        echo "[$(date)] Slopping tRNA BED ±50 nt..."
        bedtools slop \
            -i   "$REFDIR/tRNAs_fai_filtered.bed" \
            -g   {input.fai} \
            -b   50 \
            > {output.flanked_bed}
        rm -f "$REFDIR/tRNAs_fai_filtered.bed"

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
    resources:
        runtime = 60,
    conda:
        "../../envs/environment.yaml"
    shell:
        r"""
        set -euo pipefail
        exec &> {log}

        # NOTE: flanked_bed (from rule 00a) and hg38-tRNA-introns.bed both use
        # Ensembl-style chromosome names without a chr prefix — they already
        # match, so no sed conversion is needed here either.
        REFDIR=$(dirname {output.spliced_bed})
        mkdir -p "$REFDIR"

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
#
# FIX: Added `threads:` directive so Snakemake correctly reserves the cores
#      declared in params.threads from the scheduler. Previously, bowtie2-build
#      would claim cores without Snakemake's knowledge, causing over-subscription.
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
    threads: lambda wildcards: config["bowtie2"]["threads"]
    resources:
        # FIX (config tidy-up): was sge_pe="sharedmem" + a hardcoded
        # "-V -l h_vmem=4000M" literal. sge_pe collapses to 1 SGE slot
        # under the EDDIE profile's --cores 1 regardless of threads:,
        # while bowtie2-build --threads 8 genuinely uses 8 threads — same
        # class of bug fixed in trim_galore/bowtie2_pretRNA/trax_quantify.
        # Hadn't caused a visible failure yet because the pre-tRNA-only
        # reference is small enough to fit in one slot's vmem, but is
        # fixed here for consistency rather than left as a latent risk.
        runtime   = config["resources"]["build_bowtie2_pretRNA_index"]["runtime_min"],
        sge_extra = sge_extra("build_bowtie2_pretRNA_index"),
    conda:
        "../../envs/environment.yaml"
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
    resources:
        runtime = 60,
    conda:
        "../../envs/environment.yaml"
    script:
        "../scripts/build_anticodon_map.py"


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
    resources:
        runtime = 60,
    conda:
        "../../envs/environment.yaml"
    script:
        "../scripts/filter_mirbase_human.py"


# ---------------------------------------------------------------------------
# 00f-pre: Convert detailed tRNAscan-SE output to standard 9-column format
#
# The GtRNAdb tarball ships only the detailed output (15+ columns including
# HMM Score, 2' Structure Score, Isotype CM, etc.).  tRAX maketrnadb.py
# expects the classic 9-column tRNAscan-SE format:
#   Seq  tRNA#  Begin  End  Type  Anti  IntronBegin  IntronEnd  Score
# The extra columns make tRAX's parser find zero tRNA entries and abort
# with "No trna sequences".
#
# Fix: strip columns 10+ with awk, skipping the 3-line header block.
# ---------------------------------------------------------------------------
rule convert_trnascan_out:
    """
    Derive a standard 9-column tRNAscan-SE .out file from the detailed
    GtRNAdb output so that tRAX maketrnadb.py can parse it correctly.

    Two transformations applied:
      1. Strip columns 10+ (HMM Score, Isotype, etc.) — tRAX parser expects
         exactly 9 columns; extra columns cause it to return "No trna sequences".
      2. Strip the 'chr' prefix from column 1 — GtRNAdb uses UCSC-style names
         (chr1, chr2 ...) but the Ensembl primary assembly FASTA uses bare
         chromosome names (1, 2 ...).  tRAX fetches sequences from the genome
         by the name in col 1; a mismatch silently skips every entry.
         This is the same reason hg38-tRNAs_nochr.bed exists in the config.
         NOTE: convert_namemap_nochr applies the same strip to the name map so
         the tRNAscan IDs tRAX constructs (e.g. 1.trna1) continue to match.
    """
    input:
        detailed = REF["gtrndb_detailed_out"],
    output:
        standard = REF["gtrndb_standard_out"],
    log:
        f"{SCRATCH}/logs/00_convert_trnascan_out.log",
    resources:
        runtime = 60,
    shell:
        r"""
        set -euo pipefail
        exec &> {log}

        echo "[$(date)] Converting detailed tRNAscan output to standard 9-column format (chr-stripped)..."
        awk 'NR<=3{{next}} {{gsub(/^chr/,"",$1); print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}' \
            {input.detailed} > {output.standard}

        N=$(wc -l < {output.standard})
        echo "[$(date)] Done — $N tRNA entries written to {output.standard}."
        """


# ---------------------------------------------------------------------------
# 00f-pre2: Strip chr prefix from the GtRNAdb name map
#
# tRAX constructs tRNAscan IDs by joining the chromosome name from the .out
# file with the tRNA number (e.g. col1="1", col2="3" → "1.trna3").
# The raw GtRNAdb name map uses UCSC-style chr-prefixed IDs on the left
# (e.g. "chr1.trna3 → tRNA-Val-CAC-13-1").  After convert_trnascan_out strips
# "chr" from col 1 of the .out file, tRAX constructs "1.trna3" — which must
# match the left column of the name map.  This rule keeps both files in sync.
# ---------------------------------------------------------------------------
rule convert_namemap_nochr:
    """
    Strip the 'chr' prefix from tRNAscan-SE IDs in the GtRNAdb name map
    (left column) so they match the chr-stripped .out file produced by
    convert_trnascan_out.  The header line is preserved unchanged.
    """
    input:
        namemap = REF["gtrndb_name_map"],
    output:
        namemap_nochr = REF["gtrndb_name_map_nochr"],
    log:
        f"{SCRATCH}/logs/00_convert_namemap_nochr.log",
    resources:
        runtime = 60,
    shell:
        r"""
        set -euo pipefail
        exec &> {log}

        echo "[$(date)] Stripping chr prefix from name map left column..."
        awk 'NR==1{{print; next}} {{sub(/^chr/,"",$1); print $1"\t"$2}}' \
            {input.namemap} > {output.namemap_nochr}

        N=$(awk 'NR>1' {output.namemap_nochr} | wc -l)
        echo "[$(date)] Done — $N entries written to {output.namemap_nochr}."
        """


# ---------------------------------------------------------------------------
# 00f: Build TRAX database from GtRNAdb files
#      TRAX (tRAX) requires a pre-built database for tRF quantification.
#      Uses trax_env.yaml — tRAX is not in the main environment.
# ---------------------------------------------------------------------------
rule build_trax_db:
    """
    Construct the TRAX database from GtRNAdb files.

    tRAX is not a conda/pip package — it is a GitHub repository of Python
    scripts.  The trax_env.yaml conda environment provides all dependencies
    (samtools, bedtools, bowtie2, pysam …) but the scripts themselves live
    in the cloned repo at config["trax"]["script_dir"].
    All tRAX scripts are therefore invoked as:
        python {params.trax_dir}/maketrnadb.py ...
    """
    input:
        tRNA_out = REF["gtrndb_standard_out"],      # 9-column, chr-stripped (from convert_trnascan_out)
        genome   = REF["genome_fasta"],
        namemap  = REF["gtrndb_name_map_nochr"],    # chr-stripped IDs (from convert_namemap_nochr)
    output:
        db_flag  = directory(REF["trax_ref_dir"]),
    params:
        outdir   = REF["trax_ref_dir"],
        trax_dir = config["trax"]["script_dir"],
    log:
        f"{SCRATCH}/logs/00_build_trax_db.log",
    resources:
        runtime = 180,
    conda:
        "../../envs/trax_env.yaml"
    shell:
        r"""
        set -euo pipefail
        exec &> {log}
        mkdir -p {params.outdir}

        echo "[$(date)] Building TRAX database..."
        python {params.trax_dir}/maketrnadb.py \
            --trnascanfile {input.tRNA_out} \
            --genomefile   {input.genome} \
            --namemapfile  {input.namemap} \
            --databasename {params.outdir}/hsapi38

        echo "[$(date)] TRAX database built at {params.outdir}."
        """
