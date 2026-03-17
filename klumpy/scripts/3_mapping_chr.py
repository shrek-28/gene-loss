#!/usr/bin/env python3

import pandas as pd
import os
import glob
import argparse

parser = argparse.ArgumentParser(
    description="Replace chromosome numbers with RefSeq accessions"
)

parser.add_argument("-i","--input_dir",required=True)
parser.add_argument("-o","--output_dir",required=True)

parser.add_argument(
    "--map",
    nargs=2,
    action="append",
    metavar=("SPECIES","MAPFILE"),
    required=True
)

args = parser.parse_args()

INPUT_DIR = args.input_dir
OUTPUT_DIR = args.output_dir
os.makedirs(OUTPUT_DIR, exist_ok=True)

# -----------------------------
# Load mappings
# -----------------------------
species_maps = {}

for species, mapfile in args.map:

    df = pd.read_csv(mapfile, sep="\t", header=None)
    df.columns = ["RefSeq","Chromosome"]

    df["Chromosome"] = df["Chromosome"].astype(str).str.replace("SUPER_","", regex=False)

    mapping = dict(zip(df["Chromosome"], df["RefSeq"]))

    species_maps[species] = mapping

print("[INFO] Loaded mappings for:", ", ".join(species_maps.keys()))

# -----------------------------
# Process synteny tables
# -----------------------------
for file in glob.glob(os.path.join(INPUT_DIR, "*.tsv")):

    print(f"[INFO] Processing {file}")

    df = pd.read_csv(file, sep="\t")

    chr_cols = [c for c in df.columns if c.endswith("_chromosome")]

    # ensure chromosome columns are strings
    df[chr_cols] = df[chr_cols].astype(str)

    for idx, row in df.iterrows():

        species = row["Species"]

        if species not in species_maps:
            continue

        mapping = species_maps[species]

        for col in chr_cols:

            val = row[col]

            # leave untouched cases
            if (
                val.startswith("NC_") or
                val == "ABSENT" or
                val == "NA"
            ):
                continue

            if val in mapping:
                df.at[idx, col] = mapping[val]

    out_path = os.path.join(OUTPUT_DIR, os.path.basename(file))

    df.to_csv(out_path, sep="\t", index=False)

    print(f"[DONE] {out_path}")

print("[ALL DONE]")