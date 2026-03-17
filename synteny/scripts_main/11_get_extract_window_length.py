#!/usr/bin/env python3
import os
import pandas as pd
import glob
import argparse

# -----------------------------
# ARGUMENTS
# -----------------------------
parser = argparse.ArgumentParser(description="Compute extraction windows for focal gene hits.")
parser.add_argument("-b", "--best_hits_dir", required=True, help="Folder with best hit TSV files")
parser.add_argument("-n", "--neighbor_dir", required=True, help="Folder with neighbor coordinate files")
parser.add_argument("-o", "--output_dir", required=True, help="Folder to save extraction windows")
parser.add_argument("--buffer", type=int, default=20000, help="Buffer to add to windows (default 20000)")
parser.add_argument("--min_qcov", type=float, default=15, help="Minimum query coverage to include a hit (default 15)")

args = parser.parse_args()

BEST_HITS_DIR = args.best_hits_dir
NEIGHBOR_DIR = args.neighbor_dir
OUTPUT_DIR = args.output_dir
BUFFER = args.buffer
MIN_QCOV = args.min_qcov

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("[INFO] Computing extraction windows")

# -----------------------------
# LOOP OVER FOCAL GENES
# -----------------------------
for gene_dir in glob.glob(os.path.join(BEST_HITS_DIR, "*")):
    if not os.path.isdir(gene_dir):
        continue

    focal_gene = os.path.basename(gene_dir)
    print(f"[INFO] Processing {focal_gene}")

    neighbor_file = os.path.join(NEIGHBOR_DIR, f"{focal_gene}_neighbors_coordinates.tsv")
    if not os.path.exists(neighbor_file):
        print(f"[WARNING] Neighbor file missing for {focal_gene}")
        continue

    neigh_df = pd.read_csv(neighbor_file, sep="\t")
    focal_row = neigh_df[neigh_df["gene"] == focal_gene]
    if focal_row.empty:
        print(f"[WARNING] Focal gene missing in neighbors: {focal_gene}")
        continue

    # Reference spans
    S_first = neigh_df["start"].min()
    E_last = neigh_df["end"].max()
    S_focal = focal_row["start"].values[0]
    E_focal = focal_row["end"].values[0]

    # Span distances
    x = S_focal - S_first
    y = E_last - E_focal
    left_span = x + BUFFER
    right_span = y + BUFFER

    gene_out_dir = os.path.join(OUTPUT_DIR, focal_gene)
    os.makedirs(gene_out_dir, exist_ok=True)

    # -----------------------------
    # PROCESS SPECIES FILES
    # -----------------------------
    for file in glob.glob(os.path.join(gene_dir, "*_best_hit.tsv")):
        df = pd.read_csv(file, sep="\t")
        if df.empty:
            continue

        species = df["Species"].iloc[0]
        qcov = float(df["Query_Coverage"].iloc[0])

        # Skip weak/noise hits
        if qcov < MIN_QCOV:
            print(f"[SKIP] {focal_gene} {species} — low coverage")
            continue

        chrom = df["Chromosome"].iloc[0]
        hit_start = int(df["Hit_Start"].iloc[0])
        hit_end = int(df["Hit_End"].iloc[0])
        strand = df["Strand"].iloc[0]

        window_start = max(1, hit_start - left_span)
        window_end = hit_end + right_span

        out_df = pd.DataFrame([[focal_gene, species, chrom, hit_start, hit_end,
                                window_start, window_end, qcov, strand]],
                              columns=["Focal_Gene","Species","Chromosome","Hit_Start",
                                       "Hit_End","Window_Start","Window_End",
                                       "Query_Coverage","Strand"])

        out_file = os.path.join(gene_out_dir, f"{focal_gene}_{species}_extraction_window.tsv")
        out_df.to_csv(out_file, sep="\t", index=False)
        print(f"[DONE] {out_file}")

print("[ALL DONE] Extraction windows computed.")