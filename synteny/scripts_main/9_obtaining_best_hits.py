#!/usr/bin/env python3
import os
import glob
import pandas as pd
import argparse

# -----------------------------
# ARGUMENTS
# -----------------------------
parser = argparse.ArgumentParser(description="Select best focal gene hits from tblastn results.")
parser.add_argument("-t", "--tblastn_dir", required=True, help="Path to tblastn results per gene")
parser.add_argument("-n", "--neighbor_dir", required=True, help="Path to neighbor summary files")
parser.add_argument("-o", "--output_dir", required=True, help="Path to save best hit results")

args = parser.parse_args()

TBLASTN_DIR = args.tblastn_dir
NEIGHBOR_DIR = args.neighbor_dir
OUTPUT_DIR = args.output_dir

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("[INFO] Selecting best focal gene hits")

# -----------------------------
# GET STRAND FROM NEIGHBOR FILE
# -----------------------------
def get_focal_strand(gene):
    neigh_file = os.path.join(NEIGHBOR_DIR, f"{gene}_neighbors_coordinates.tsv")
    if not os.path.exists(neigh_file):
        return "+"
    df = pd.read_csv(neigh_file, sep="\t")
    row = df[df["gene"] == gene]
    if row.empty:
        return "+"
    return row.iloc[0]["strand"]

# -----------------------------
# LOOP OVER GENES
# -----------------------------
for gene_dir in glob.glob(os.path.join(TBLASTN_DIR, "*")):
    if not os.path.isdir(gene_dir):
        continue

    focal_gene = os.path.basename(gene_dir)
    print(f"\n[INFO] Processing {focal_gene}")

    strand = get_focal_strand(focal_gene)

    gene_out_dir = os.path.join(OUTPUT_DIR, focal_gene)
    os.makedirs(gene_out_dir, exist_ok=True)

    # -----------------------------
    # LOOP OVER SPECIES
    # -----------------------------
    for tsv_file in glob.glob(os.path.join(gene_dir, "*.tsv")):
        if tsv_file.endswith(".best.tsv"):
            continue

        species = os.path.basename(tsv_file).replace(f"{focal_gene}_vs_", "").replace(".tsv", "")
        if os.path.getsize(tsv_file) == 0:
            continue

        cols = ["qseqid","sseqid","sstart","send","evalue","bitscore","qcovs"]
        df = pd.read_csv(tsv_file, sep="\t", header=None, names=cols)
        if df.empty:
            continue

        for col in ["bitscore","qcovs","sstart","send"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        df = df.dropna(subset=["bitscore","qcovs","sstart","send"])
        if df.empty:
            continue

        # -----------------------------
        # KEEP CHROMOSOMES ONLY
        # -----------------------------
        # Remove 'ref|' prefix and trailing '|'
        df["sseqid"] = df["sseqid"].str.replace(r"ref\|(.+?)\|", r"\1", regex=True)

        # Then keep only NC chromosomes
        df = df[df["sseqid"].str.startswith("NC_")]
        if df.empty:
            continue

        df["hit_start"] = df[["sstart","send"]].min(axis=1)
        df["hit_end"]   = df[["sstart","send"]].max(axis=1)

        # -----------------------------
        # GROUP BY CHROMOSOME
        # -----------------------------
        chrom_grp = df.groupby("sseqid").agg({"qcovs":"max","bitscore":"max"}).reset_index()
        chrom_grp = chrom_grp.sort_values(["qcovs","bitscore"], ascending=False)
        best_chr = chrom_grp.iloc[0]["sseqid"]

        # -----------------------------
        # BEST HIT WITHIN CHROMOSOME
        # -----------------------------
        sub = df[df["sseqid"] == best_chr]
        sub = sub.sort_values(["qcovs","bitscore"], ascending=False)
        best_hit = sub.iloc[0]

        # -----------------------------
        # OUTPUT
        # -----------------------------
        out_df = pd.DataFrame([{
            "Focal_Gene": focal_gene,
            "Species": species,
            "Chromosome": best_chr,
            "Hit_Start": int(best_hit["hit_start"]),
            "Hit_End": int(best_hit["hit_end"]),
            "Query_Coverage": float(best_hit["qcovs"]),
            "Bitscore": float(best_hit["bitscore"]),
            "Strand": strand
        }])

        out_file = os.path.join(gene_out_dir, f"{focal_gene}_{species}_best_hit.tsv")
        out_df.to_csv(out_file, sep="\t", index=False)
        print(f"[OK] {species}")

print("\n[DONE] Best focal gene hits generated.")