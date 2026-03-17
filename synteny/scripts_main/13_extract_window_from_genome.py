#!/usr/bin/env python3
import os
import pandas as pd
import subprocess
import glob
import argparse

# -----------------------------
# ARGUMENT PARSER
# -----------------------------
parser = argparse.ArgumentParser(
    description="Extract genomic regions using samtools faidx based on extraction window TSV files."
)

parser.add_argument("--window_root", required=True, help="Directory containing extraction window folders")
parser.add_argument("--genome_dir", required=True, help="Directory containing cleaned genome FASTA files")
parser.add_argument("--output_root", required=True, help="Directory where extracted regions will be written")

args = parser.parse_args()

WINDOW_ROOT = args.window_root
GENOME_DIR = args.genome_dir
OUTPUT_ROOT = args.output_root

os.makedirs(OUTPUT_ROOT, exist_ok=True)

# -----------------------------
# LOOP OVER FOCAL GENES
# -----------------------------
for gene_dir in glob.glob(os.path.join(WINDOW_ROOT, "*")):

    if not os.path.isdir(gene_dir):
        continue

    focal_gene = os.path.basename(gene_dir)
    print(f"\n[INFO] Processing focal gene: {focal_gene}")

    gene_out_dir = os.path.join(OUTPUT_ROOT, focal_gene)
    os.makedirs(gene_out_dir, exist_ok=True)

    # -----------------------------
    # LOOP OVER SPECIES FILES
    # -----------------------------
    for tsv_file in glob.glob(os.path.join(gene_dir, "*_extraction_window.tsv")):

        df = pd.read_csv(tsv_file, sep="\t")

        if df.empty:
            continue

        row = df.iloc[0]

        species = row["Species"]
        chrom = row["Chromosome"]
        start = int(row["Window_Start"])
        end = int(row["Window_End"])

        genome_file = os.path.join(
            GENOME_DIR,
            f"{species}.fna"
        )

        if not os.path.exists(genome_file):
            print(f"[WARNING] Genome missing: {species}")
            continue

        species_out_dir = os.path.join(gene_out_dir, species)
        os.makedirs(species_out_dir, exist_ok=True)

        out_fasta = os.path.join(
            species_out_dir,
            f"{focal_gene}_{species}_region.fna"
        )

        region = f"{chrom}:{start}-{end}"

        try:
            with open(out_fasta, "w") as out_handle:
                subprocess.run(
                    ["samtools", "faidx", genome_file, region],
                    check=True,
                    stdout=out_handle
                )

            print(f"[OK] {focal_gene} | {species}")

        except subprocess.CalledProcessError:
            print(f"[FAILED] {focal_gene} | {species}")

print("\n[DONE] All genomic regions extracted.")