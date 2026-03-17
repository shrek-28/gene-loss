#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
import glob
import argparse
from matplotlib.lines import Line2D

# --------------------------------------------------
# Argument parser
# --------------------------------------------------
parser = argparse.ArgumentParser(
    description="Plot gene presence/absence across species using BLAST results."
)

parser.add_argument(
    "--neigh_coord_dir",
    required=True,
    help="Directory containing *_neighbors_coordinates.tsv files"
)

parser.add_argument(
    "--coord_root",
    required=True,
    help="Root directory containing neighbouring_genes_genomic_coordinates"
)

parser.add_argument(
    "--species_list",
    required=True,
    help="Text file containing species names (one per line)"
)

parser.add_argument(
    "--out_dir",
    required=True,
    help="Output directory for plots"
)

args = parser.parse_args()

NEIGH_COORD_DIR = args.neigh_coord_dir
COORD_ROOT = args.coord_root
SPECIES_LIST = args.species_list
OUT_DIR = args.out_dir

os.makedirs(OUT_DIR, exist_ok=True)

# --------------------------------------------------
# Load species list
# --------------------------------------------------
with open(SPECIES_LIST) as f:
    species_list = [x.strip() for x in f if x.strip()]

# --------------------------------------------------
# Iterate focal genes
# --------------------------------------------------
for coord_file in glob.glob(os.path.join(
        NEIGH_COORD_DIR,
        "*_neighbors_coordinates.tsv")):

    focal_gene = os.path.basename(coord_file).replace(
        "_neighbors_coordinates.tsv", "")

    print(f"[INFO] Plotting {focal_gene}")

    df = pd.read_csv(coord_file, sep="\t")

    genes = df.sort_values("start")["gene"].tolist()

    gene_dirs = []
    for g in genes:
        if g == focal_gene:
            gene_dirs.append(f"{g}_focalgene")
        else:
            gene_dirs.append(g)

    focal_root = os.path.join(COORD_ROOT, focal_gene)
    if not os.path.exists(focal_root):
        print("  [SKIP] Missing coordinate folder")
        continue

    # --------------------------------------------------
    # Plot
    # --------------------------------------------------
    plt.figure(figsize=(18, 11))

    for y, species in enumerate(species_list):

        for x, gene_dir in enumerate(gene_dirs):

            gene_path = os.path.join(focal_root, gene_dir)

            if not os.path.exists(gene_path):
                continue

            matches = glob.glob(
                os.path.join(gene_path, f"*vs_{species}.best.tsv")
            )

            if len(matches) == 0:
                continue

            best_file = matches[0]

            if os.path.getsize(best_file) == 0:
                plt.scatter(x, y,
                            facecolors='none',
                            edgecolors='black',
                            s=90)
                continue

            df_hit = pd.read_csv(best_file, sep="\t")

            if df_hit.empty:
                plt.scatter(x, y,
                            facecolors='none',
                            edgecolors='black',
                            s=90)
                continue

            strand = df_hit.get("Strand", ["+"])[0]

            if strand == "+":
                color = "pink"
            else:
                color = "blue"

            if gene_dir.endswith("_focalgene"):
                color = "red" if strand == "+" else "navy"

            plt.scatter(x, y, color=color, s=110)

    # --------------------------------------------------
    # Axis formatting
    # --------------------------------------------------
    plt.yticks(range(len(species_list)), species_list,
               fontsize=12, fontweight='bold')

    plt.xticks(range(len(genes)), genes,
               rotation=90,
               fontsize=12,
               fontweight='bold')

    plt.xlabel("Gene order (Homo sapiens reference)",
               fontsize=14,
               fontweight='bold')

    plt.ylabel("Species",
               fontsize=14,
               fontweight='bold')

    plt.title(f"Syntenic gene presence — {focal_gene}",
              fontsize=18,
              fontweight='bold')

    # --------------------------------------------------
    # Legend
    # --------------------------------------------------
    legend_items = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor='pink', markersize=8,
               label='Gene (+ strand)'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor='blue', markersize=8,
               label='Gene (− strand)'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor='red', markersize=8,
               label='Focal gene (+ strand)'),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor='navy', markersize=8,
               label='Focal gene (− strand)'),
        Line2D([0], [0], marker='o', color='black',
               markerfacecolor='none', markersize=8,
               label='No BLAST hit'),
        Line2D([], [], linestyle='None',
               label='No circle = gene absent')
    ]

    plt.legend(handles=legend_items,
               bbox_to_anchor=(1.02, 1),
               loc="upper left",
               prop={'weight': 'bold', 'size': 12})

    # --------------------------------------------------
    # Save
    # --------------------------------------------------
    out_file = os.path.join(
        OUT_DIR,
        f"{focal_gene}_gene_presence.png"
    )

    plt.tight_layout()
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"  [SAVED] {out_file}")

print("\n[DONE] All plots generated.")