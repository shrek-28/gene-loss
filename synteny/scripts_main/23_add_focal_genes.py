#!/usr/bin/env python3
import os
import shutil
import glob
import argparse

parser = argparse.ArgumentParser(description="Copy focal gene best hits into neighbour structure.")
parser.add_argument("-i", "--input", required=True, help="Input directory containing focal gene best hits")
parser.add_argument("-o", "--output", required=True, help="Output root directory for copied files")
args = parser.parse_args()

FOCAL_INPUT = args.input
OUTPUT_ROOT = args.output

print("[INFO] Copying focal genes into neighbour structure")

for focal_dir in glob.glob(os.path.join(FOCAL_INPUT, "*")):

    if not os.path.isdir(focal_dir):
        continue

    focal_gene = os.path.basename(focal_dir)

    target_root = os.path.join(OUTPUT_ROOT, focal_gene)
    os.makedirs(target_root, exist_ok=True)

    # create focal gene folder with suffix
    focal_out_dir = os.path.join(target_root, f"{focal_gene}_focalgene")
    os.makedirs(focal_out_dir, exist_ok=True)

    for file in glob.glob(os.path.join(focal_dir, "*_best_hit.tsv")):

        fname = os.path.basename(file)

        species = fname.replace(f"{focal_gene}_", "").replace("_best_hit.tsv", "")

        new_name = f"{focal_gene}_{focal_gene}_vs_{species}.best.tsv"
        out_file = os.path.join(focal_out_dir, new_name)

        shutil.copy(file, out_file)

        print(f"[OK] {focal_gene} | {species}")

print("[DONE] Focal genes copied")