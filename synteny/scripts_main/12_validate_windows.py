#!/usr/bin/env python3
import os
import pandas as pd
import argparse

# -----------------------------
# ARGUMENTS
# -----------------------------
parser = argparse.ArgumentParser(description="Validate extraction windows against neighbor reference spans.")
parser.add_argument("-n", "--neighbor_dir", required=True, help="Folder with neighbor coordinate files")
parser.add_argument("-w", "--window_dir", required=True, help="Folder with extraction windows")
parser.add_argument("-s", "--species", required=True, help="Species name to check in the window filenames (e.g., 'Homo_sapiens')")

args = parser.parse_args()

NEIGHBOR_DIR = args.neighbor_dir
WINDOW_DIR = args.window_dir
SPECIES = args.species

print(f"\nChecking {SPECIES} extraction windows\n")

# -----------------------------
# LOOP OVER FOCAL GENES
# -----------------------------
for gene in os.listdir(WINDOW_DIR):
    gene_dir = os.path.join(WINDOW_DIR, gene)
    if not os.path.isdir(gene_dir):
        continue

    # ---------- neighbor reference ----------
    neigh_file = os.path.join(NEIGHBOR_DIR, f"{gene}_neighbors_coordinates.tsv")
    if not os.path.exists(neigh_file):
        print(f"[SKIP] No neighbor file for {gene}")
        continue

    neigh_df = pd.read_csv(neigh_file, sep="\t")
    ref_start = neigh_df["start"].min()
    ref_end = neigh_df["end"].max()
    ref_span = ref_end - ref_start

    # ---------- species window ----------
    win_file = os.path.join(gene_dir, f"{gene}_{SPECIES}_extraction_window.tsv")
    if not os.path.exists(win_file):
        print(f"[SKIP] No {SPECIES} window for {gene}")
        continue

    win_df = pd.read_csv(win_file, sep="\t")
    win_start = win_df["Window_Start"].iloc[0]
    win_end = win_df["Window_End"].iloc[0]
    win_span = win_end - win_start

    # ---------- comparison ----------
    if win_span < ref_span:
        print(f"[PROBLEM] {gene}")
        print(f"  Reference span : {ref_span}")
        print(f"  Extracted span : {win_span}")
    else:
        print(f"[OK] {gene}")
        print(f"  Reference span : {ref_span}")
        print(f"  Extracted span : {win_span}")

print("\nValidation complete\n")