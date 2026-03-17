#!/usr/bin/env python3
import os
import pandas as pd
import glob
import argparse

# -----------------------------
# Argument parsing
# -----------------------------
parser = argparse.ArgumentParser(description="Build clean synteny tables from neighbour genomic coordinates.")
parser.add_argument("-c", "--coord_dir", required=True, help="Directory containing neighbour summary files")
parser.add_argument("-g", "--gene_coord_root", required=True, help="Root directory for neighbouring genes genomic coordinates")
parser.add_argument("-s", "--species_list", required=True, help="Path to species list file")
parser.add_argument("-o", "--out_root", required=True, help="Output directory for synteny tables")
args = parser.parse_args()

COORD_DIR = args.coord_dir
GENE_COORD_ROOT = args.gene_coord_root
SPECIES_LIST = args.species_list
OUT_ROOT = args.out_root
os.makedirs(OUT_ROOT, exist_ok=True)

print("[INFO] Building clean synteny tables")

# -----------------------------
# Load species list
# -----------------------------
with open(SPECIES_LIST) as f:
    species_list = [l.strip() for l in f if l.strip()]

# -----------------------------
# Loop focal genes
# -----------------------------
for coord_file in os.listdir(COORD_DIR):

    if not coord_file.endswith("_neighbors_coordinates.tsv"):
        continue

    focal_gene = coord_file.replace("_neighbors_coordinates.tsv", "")
    print(f"[INFO] Processing {focal_gene}")

    coord_path = os.path.join(COORD_DIR, coord_file)

    ref_df = pd.read_csv(coord_path, sep="\t")
    ref_df = ref_df.sort_values("start")

    gene_order = ref_df["gene"].tolist()

    rows = []
    gene_dir = os.path.join(GENE_COORD_ROOT, focal_gene)

    for species in species_list:

        row = {"Species": species}

        for gene in gene_order:

            # focal gene stored differently
            if gene == focal_gene:
                gene_subdir = os.path.join(
                    gene_dir,
                    f"{gene}_focalgene"
                )
            else:
                gene_subdir = os.path.join(gene_dir, gene)

            if not os.path.exists(gene_subdir):
                # row[gene] = "ABSENT"
                row[f"{gene}_start"] = "ABSENT"
                row[f"{gene}_end"] = "ABSENT"
                row[f"{gene}_strand"] = "ABSENT"
                row[f"{gene}_chromosome"] = "ABSENT"
                continue

            # strict match for species
            pattern = os.path.join(
                gene_subdir,
                f"*_{species}.best.tsv"
            )

            matches = glob.glob(pattern)

            if len(matches) == 0:
                row[f"{gene}_start"] = "ABSENT"
                row[f"{gene}_end"] = "ABSENT"
                row[f"{gene}_strand"] = "ABSENT"
                row[f"{gene}_chromosome"] = "ABSENT"
                continue

            file_path = matches[0]

            if os.path.getsize(file_path) == 0:
                # row[gene] = "NA"
                row[f"{gene}_start"] = "NA"
                row[f"{gene}_end"] = "NA"
                row[f"{gene}_strand"] = "NA"
                row[f"{gene}_chromosome"] = "NA"
                continue

            df = pd.read_csv(file_path, sep="\t")

            if df.empty:
                # row[gene] = "NA"
                row[f"{gene}_start"] = "NA"
                row[f"{gene}_end"] = "NA"
                row[f"{gene}_strand"] = "NA"
                row[f"{gene}_chromosome"] = "NA"
                continue

            # prefer genomic coordinates
            # if "Genomic_Start" in df.columns:
            #     row[gene] = int(df.iloc[0]["Genomic_Start"])
            # else:
            #     row[gene] = int(df.iloc[0]["Hit_Start"])
            if "Genomic_Start" in df.columns:
                row[f"{gene}_start"] = int(df.iloc[0]["Genomic_Start"])
                row[f"{gene}_end"] = int(df.iloc[0]["Genomic_End"])
                row[f"{gene}_strand"] = df.iloc[0]["Strand"]
                row[f"{gene}_chromosome"] = df.iloc[0]["Chromosome"]
            else:
                row[f"{gene}_start"] = int(df.iloc[0]["Hit_Start"])
                row[f"{gene}_end"] = int(df.iloc[0]["Hit_End"])
                row[f"{gene}_strand"] = df.iloc[0]["Strand"]
                row[f"{gene}_chromosome"] = df.iloc[0]["Chromosome"]

        rows.append(row)

    out_df = pd.DataFrame(rows)

    cols = ["Species"]

    for gene in gene_order:
        cols.extend([
            f"{gene}_start",
            f"{gene}_end",
            f"{gene}_strand",
            f"{gene}_chromosome"
        ])
    out_df = out_df[cols]

    out_file = os.path.join(
        OUT_ROOT,
        f"{focal_gene}_synteny_table.tsv"
    )

    out_df.to_csv(out_file, sep="\t", index=False)

    print(f"[DONE] {out_file}")

print("[ALL DONE] Clean synteny tables created.")