import logging

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO,
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

n_mat = filter_fasta(snakemake.input.mature,  snakemake.output.mature_hsa)
n_hp  = filter_fasta(snakemake.input.hairpin, snakemake.output.hairpin_hsa)
logging.info(f"Kept {n_mat} mature hsa miRNAs, {n_hp} hairpin hsa miRNAs.")
