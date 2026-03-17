#!/usr/bin/env python3
import os
import pandas as pd
import glob
import argparse

# -----------------------------
# ARGUMENTS
# -----------------------------
parser = argparse.ArgumentParser(description="Add gene presence status based on query coverage.")
parser.add_argument("-i", "--input_dir", required=True, help="Path to the folder containing best hit TSV files")
parser.add_argument("-o", "--output_dir", required=True, help="Path to save status-annotated TSV files")
args = parser.parse_args()

BEST_HITS_DIR = args.input_dir
OUTPUT_DIR = args.output_dir

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("[INFO] Adding gene presence status based on query coverage")

# -----------------------------
# Status classification
# -----------------------------
def classify_cov(cov):
    if cov >= 50:
        return "Strong_gene_presence"
    elif cov >= 30:
        return "Partial_degraded_gene"
    elif cov >= 15:
        return "Weak_evidence_warning"
    else:
        return "Likely_noise"

# -----------------------------
# Process files
# -----------------------------
for gene_dir in glob.glob(os.path.join(BEST_HITS_DIR, "*")):
    if not os.path.isdir(gene_dir):
        continue

    gene = os.path.basename(gene_dir)
    out_gene_dir = os.path.join(OUTPUT_DIR, gene)
    os.makedirs(out_gene_dir, exist_ok=True)

    for file in glob.glob(os.path.join(gene_dir, "*.tsv")):
        df = pd.read_csv(file, sep="\t")

        # ensure coverage column exists
        if "Query_Coverage" not in df.columns:
            print(f"[WARNING] Query_Coverage missing in {file}")
            continue

        df["Status"] = df["Query_Coverage"].apply(classify_cov)

        out_file = os.path.join(out_gene_dir, os.path.basename(file))
        df.to_csv(out_file, sep="\t", index=False)
        print(f"[DONE] {out_file}")

print("[ALL DONE] Status assignment completed.")