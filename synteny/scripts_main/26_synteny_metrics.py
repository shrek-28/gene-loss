#!/usr/bin/env python3
import os
import pandas as pd
import glob
import argparse

# -----------------------------------
# Argument parser
# -----------------------------------
parser = argparse.ArgumentParser(
    description="Compute synteny conservation metrics from synteny tables."
)

parser.add_argument(
    "--table_dir",
    required=True,
    help="Directory containing *_synteny_table.tsv files"
)

parser.add_argument(
    "--ref_dir",
    required=True,
    help="Directory containing *_neighbors_coordinates.tsv reference files"
)

parser.add_argument(
    "--out_dir",
    required=True,
    help="Directory where metrics will be written"
)

args = parser.parse_args()

TABLE_DIR = args.table_dir
REF_DIR = args.ref_dir
OUT_DIR = args.out_dir

os.makedirs(OUT_DIR, exist_ok=True)

print("[INFO] Computing synteny conservation metrics")

# -----------------------------------
# PROCESS EACH FOCAL GENE
# -----------------------------------
for table_file in glob.glob(os.path.join(TABLE_DIR, "*_synteny_table.tsv")):

    focal_gene = os.path.basename(table_file).replace(
        "_synteny_table.tsv", ""
    )

    print(f"[INFO] Processing {focal_gene}")

    df = pd.read_csv(table_file, sep="\t")
    genes = df.columns[1:].tolist()

    if len(genes) == 0:
        print(f"[WARNING] No gene columns found in {table_file}. Skipping.")
        continue
    # -----------------------------------
    # Load reference gene order
    # -----------------------------------
    ref_file = os.path.join(
        REF_DIR,
        f"{focal_gene}_neighbors_coordinates.tsv"
    )

    ref_df = pd.read_csv(ref_file, sep="\t")
    ref_df = ref_df.sort_values("start")

    ref_genes = ref_df["gene"].tolist()

    ref_pairs = [
        (ref_genes[i], ref_genes[i + 1])
        for i in range(len(ref_genes) - 1)
    ]

    results = []

    # -----------------------------------
    # Score each species
    # -----------------------------------
    for _, row in df.iterrows():

        species = row["Species"]

        # ---------- Presence ----------
        present_genes = []
        for g in genes:
            val = row[g]
            if val not in ["ABSENT", "NA"]:
                present_genes.append(g)

        presence_score = len(present_genes) / len(genes)

        # ---------- Coordinates ----------
        coords = {}
        for g in present_genes:
            coords[g] = float(row[g])

        species_order = sorted(coords.items(), key=lambda x: x[1])
        species_order = [x[0] for x in species_order]

        # ---------- Adjacency ----------
        valid_pairs = 0
        correct_pairs = 0

        for g1, g2 in ref_pairs:

            if g1 in coords and g2 in coords:
                valid_pairs += 1

                idx1 = species_order.index(g1)
                idx2 = species_order.index(g2)

                if abs(idx1 - idx2) == 1:
                    correct_pairs += 1

        if valid_pairs > 0:
            order_score = correct_pairs / valid_pairs
        else:
            order_score = 0

        # ---------- Final ----------
        final_score = 0.5 * presence_score + 0.5 * order_score

        results.append([
            species,
            round(presence_score, 4),
            round(order_score, 4),
            round(final_score, 4)
        ])

    out_df = pd.DataFrame(
        results,
        columns=[
            "Species",
            "Presence_Score",
            "Order_Score",
            "Synteny_Score"
        ]
    )

    out_file = os.path.join(
        OUT_DIR,
        f"{focal_gene}_metrics.tsv"
    )

    out_df.to_csv(out_file, sep="\t", index=False)

    print(f"[OK] Saved {out_file}")

print("\n[DONE] Metrics generated per gene")