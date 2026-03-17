#!/usr/bin/env python3
import os
import pandas as pd
import glob
import argparse

# ------------------------------------------------
# ARGUMENT PARSER
# ------------------------------------------------
parser = argparse.ArgumentParser(
    description="Select best BLAST hits for non-protein-coding genes."
)

parser.add_argument(
    "--root_dir",
    required=True,
    help="Root directory containing blast_non_proteincoding results"
)

args = parser.parse_args()
ROOT_DIR = args.root_dir

print("[INFO] Selecting best hits for non-protein coding genes")

# ------------------------------------------------
# Status classification
# ------------------------------------------------
def classify_status(qcov):
    if qcov >= 50:
        return "STRONG"
    elif qcov >= 30:
        return "PARTIAL"
    elif qcov >= 15:
        return "WEAK"
    else:
        return "NOISE"

# ------------------------------------------------
# Loop structure
# ------------------------------------------------
for focal_gene_dir in glob.glob(os.path.join(ROOT_DIR, "*")):

    if not os.path.isdir(focal_gene_dir):
        continue

    focal_gene = os.path.basename(focal_gene_dir)

    for gene_dir in glob.glob(os.path.join(focal_gene_dir, "*")):

        if not os.path.isdir(gene_dir):
            continue

        gene_id = os.path.basename(gene_dir)

        for species_dir in glob.glob(os.path.join(gene_dir, "*")):

            if not os.path.isdir(species_dir):
                continue

            species = os.path.basename(species_dir)

            for tsv_file in glob.glob(os.path.join(species_dir, "*.blastn.tsv")):

                if os.path.getsize(tsv_file) == 0:
                    continue

                # Load BLAST output
                df = pd.read_csv(tsv_file, sep="\t", header=None)

                df.columns = [
                    "qseqid","sseqid","pident","length",
                    "mismatch","gapopen","qstart","qend",
                    "sstart","send","evalue","bitscore","qcovs"
                ]

                # Ensure numeric columns
                for col in ["pident","bitscore","qcovs","sstart","send"]:
                    df[col] = pd.to_numeric(df[col], errors="coerce")

                df = df.dropna(subset=["qcovs","bitscore"])
                if df.empty:
                    continue

                # Normalize coordinates
                df["hit_start"] = df[["sstart","send"]].min(axis=1)
                df["hit_end"] = df[["sstart","send"]].max(axis=1)

                # Strand orientation
                df["strand"] = df.apply(
                    lambda r: "-" if r["send"] < r["sstart"] else "+",
                    axis=1
                )

                # Select best hit
                df = df.sort_values(
                    ["qcovs","bitscore","pident"],
                    ascending=False
                )

                best = df.iloc[0]
                status = classify_status(float(best["qcovs"]))

                best_df = pd.DataFrame([{
                    "Focal_Gene": focal_gene,
                    "Gene_ID": gene_id,
                    "Species": species,
                    "Chromosome": best["sseqid"],
                    "Hit_Start": int(best["hit_start"]),
                    "Hit_End": int(best["hit_end"]),
                    "Strand": best["strand"],
                    "Query_Coverage": float(best["qcovs"]),
                    "Bitscore": float(best["bitscore"]),
                    "Percent_Identity": float(best["pident"]),
                    "Status": status
                }])

                out_file = tsv_file.replace(".blastn.tsv", ".best.tsv")
                best_df.to_csv(out_file, sep="\t", index=False)

                print(f"[OK] {focal_gene} | {gene_id} | {species}")

print("[DONE] Best hits selected for non-protein coding genes")