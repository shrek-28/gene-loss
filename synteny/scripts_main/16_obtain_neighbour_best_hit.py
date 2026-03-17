#!/usr/bin/env python3
import os
import pandas as pd
import glob
import argparse

# -----------------------------
# ARGUMENT PARSER
# -----------------------------
parser = argparse.ArgumentParser(
    description="Select best tblastn hits and assign strength status."
)

parser.add_argument(
    "--root_dir",
    required=True,
    help="Directory containing tblastn results for neighbouring genes"
)

args = parser.parse_args()
ROOT_DIR = args.root_dir

print("[INFO] Selecting best hits and assigning status")

# -----------------------------
# Status classification
# -----------------------------
def classify_status(qcov):
    if qcov >= 50:
        return "STRONG"
    elif qcov >= 30:
        return "PARTIAL"
    elif qcov >= 15:
        return "WEAK"
    else:
        return "NOISE"

# -----------------------------
# LOOP STRUCTURE
# -----------------------------
for focal_gene_dir in glob.glob(os.path.join(ROOT_DIR, "*")):

    if not os.path.isdir(focal_gene_dir):
        continue

    focal_gene = os.path.basename(focal_gene_dir)

    for neigh_gene_dir in glob.glob(os.path.join(focal_gene_dir, "*")):

        if not os.path.isdir(neigh_gene_dir):
            continue

        neigh_gene = os.path.basename(neigh_gene_dir)

        for tsv_file in glob.glob(os.path.join(neigh_gene_dir, "*.tblastn.tsv")):

            if os.path.getsize(tsv_file) == 0:
                continue

            species = os.path.basename(tsv_file).split("_vs_")[-1].replace(".tblastn.tsv", "")

            df = pd.read_csv(tsv_file, sep="\t", header=None)

            df.columns = [
                "qseqid","sseqid","pident","length",
                "mismatch","gapopen","qstart","qend",
                "sstart","send","evalue","bitscore","qcovs"
            ]

            # normalize hit coords
            df["hit_start"] = df[["sstart","send"]].min(axis=1)
            df["hit_end"] = df[["sstart","send"]].max(axis=1)

            # strand from BLAST orientation
            df["strand"] = df.apply(
                lambda r: "-" if r["send"] < r["sstart"] else "+",
                axis=1
            )

            # choose best hit
            df = df.sort_values(
                ["qcovs", "bitscore", "pident"],
                ascending=False
            )

            best = df.iloc[0]

            status = classify_status(best["qcovs"])

            best_df = pd.DataFrame([{
                "Focal_Gene": focal_gene,
                "Neighbour_Gene": neigh_gene,
                "Species": species,
                "Chromosome": best["sseqid"],
                "Hit_Start": int(best["hit_start"]),
                "Hit_End": int(best["hit_end"]),
                "Strand": best["strand"],
                "Query_Coverage": best["qcovs"],
                "Bitscore": best["bitscore"],
                "Percent_Identity": best["pident"],
                "Status": status
            }])

            out_file = tsv_file.replace(".tblastn.tsv", ".best.tsv")
            best_df.to_csv(out_file, sep="\t", index=False)

            print(f"[OK] {focal_gene} | {neigh_gene} | {species}")

print("[DONE] Best hits extracted with status")