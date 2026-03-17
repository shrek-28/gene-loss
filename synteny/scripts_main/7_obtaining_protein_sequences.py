#!/usr/bin/env python3
import os
import argparse
from pathlib import Path
from Bio import SeqIO

# -------------------------------
# Argument parser
# -------------------------------
parser = argparse.ArgumentParser(
    description="Extract protein FASTA files per focal gene neighborhood."
)

parser.add_argument("--tsv_dir", required=True,
                    help="Directory containing *_neighbour_summary.tsv files")
parser.add_argument("--proteins", required=True,
                    help="Path to cleaned protein FASTA (.faa)")
parser.add_argument("--output_dir", required=True,
                    help="Directory to write per-gene FASTA outputs")

args = parser.parse_args()

TSV_DIR = Path(args.tsv_dir)
PROTEIN_FASTA = Path(args.proteins)
OUTPUT_DIR = Path(args.output_dir)
OUTPUT_DIR.mkdir(exist_ok=True)

# -------------------------------
# Load protein FASTA
# -------------------------------
print("[INFO] Loading protein FASTA...")
protein_dict = {rec.id: rec for rec in SeqIO.parse(PROTEIN_FASTA, "fasta")}
print(f"[INFO] {len(protein_dict)} proteins loaded.")

# -------------------------------
# Process each focal gene TSV
# -------------------------------
for tsv_file in TSV_DIR.glob("*_neighbour_summary.tsv"):
    focal_gene = tsv_file.stem.replace("_neighbour_summary", "")
    print(f"[INFO] Processing focal gene: {focal_gene}")

    gene_dir = OUTPUT_DIR / focal_gene
    gene_dir.mkdir(exist_ok=True)

    with open(tsv_file) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                continue

            gene_id, gene_type, longest_protein = parts

            if gene_type == "protein_coding" and longest_protein != "NA":
                if longest_protein in protein_dict:
                    out_file = gene_dir / f"{gene_id}.faa"
                    SeqIO.write(protein_dict[longest_protein], out_file, "fasta")
                    print(f"[OK] {gene_id} -> {longest_protein}")
                else:
                    print(f"[WARNING] Protein {longest_protein} not found for {gene_id}")
            else:
                out_file = gene_dir / f"{gene_id}_status.txt"
                with open(out_file, "w") as o:
                    o.write(f"{gene_type}\n")
                print(f"[INFO] {gene_id} is {gene_type}")

print("[DONE] All focal genes processed. Results in", OUTPUT_DIR)