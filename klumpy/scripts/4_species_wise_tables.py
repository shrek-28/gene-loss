#!/usr/bin/env python3

import pandas as pd
import os
import glob
import argparse

parser = argparse.ArgumentParser(
    description="Generate species boundary tables from synteny tables"
)

parser.add_argument(
    "-i",
    "--input_dir",
    required=True,
    help="Directory containing *_synteny_table.tsv files"
)

parser.add_argument(
    "-o",
    "--output_dir",
    required=True,
    help="Output directory"
)

args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = args.output_dir

os.makedirs(OUTPUT_DIR, exist_ok=True)

species_results = {}

# -----------------------------
# Process each synteny table
# -----------------------------
for file in glob.glob(os.path.join(INPUT_DIR, "*_synteny_table.tsv")):

    focal_gene = os.path.basename(file).replace("_synteny_table.tsv", "")

    df = pd.read_csv(file, sep="\t")

    start_cols = [f"{focal_gene}_start"] if f"{focal_gene}_start" in df.columns else []
    end_cols   = [f"{focal_gene}_end"]   if f"{focal_gene}_end"   in df.columns else []

    for _, row in df.iterrows():

        species = row["Species"]

        # collect numeric start values
        starts = []
        for col in start_cols:
            val = row[col]
            if val != "ABSENT" and val != "NA":
                starts.append(int(val))

        # collect numeric end values
        ends = []
        for col in end_cols:
            val = row[col]
            if val != "ABSENT" and val != "NA":
                ends.append(int(val))

        if not starts or not ends:
            continue

        leftbound = min(starts)
        rightbound = max(ends)

        chr_col = f"{focal_gene}_chromosome"

        chromosome = row[chr_col]

        entry = {
            "focal_gene": focal_gene,
            "chromosome": chromosome,
            "leftbound": leftbound,
            "rightbound": rightbound
        }

        if species not in species_results:
            species_results[species] = []

        species_results[species].append(entry)

# -----------------------------
# Write species tables
# -----------------------------
for species, rows in species_results.items():

    out_df = pd.DataFrame(rows)

    out_file = os.path.join(OUTPUT_DIR, f"{species}.tsv")

    out_df.to_csv(out_file, sep="\t", index=False)

    print(f"[DONE] {out_file}")

print("[ALL DONE]")