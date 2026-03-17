#!/usr/bin/env python3

import argparse
import pandas as pd
import subprocess
import os
import glob
import pysam
import shutil

parser = argparse.ArgumentParser(description="Generate alignment plots for all genes")

parser.add_argument("--bam", required=True, help="Sorted BAM file")
parser.add_argument("--tsv", required=True, help="TSV with gene regions")
parser.add_argument("--outdir", required=True, help="Directory for PDF plots")
parser.add_argument("--table_dir", required=True, help="Directory for gene tables")
parser.add_argument("--summary", required=True, help="Output flanking summary TSV")
parser.add_argument("--annotation", required=False, help="add annotation file for genes")

args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
os.makedirs(args.table_dir, exist_ok=True)
os.makedirs(os.path.dirname(args.summary), exist_ok=True)

df = pd.read_csv(args.tsv, sep="\t")
bam = pysam.AlignmentFile(args.bam, "rb")

flank_sizes = [500000, 250000, 100000, 50000]

summary = []

for _, row in df.iterrows():

    gene = row["focal_gene"]
    chrom = row["chromosome"]
    gene_start = int(row["leftbound"])
    gene_end = int(row["rightbound"])

    chrom_len = bam.get_reference_length(chrom)

    # gene-specific directories
    gene_plot_dir = os.path.join(args.outdir, gene)
    gene_table_dir = os.path.join(args.table_dir, gene)

    os.makedirs(gene_plot_dir, exist_ok=True)
    os.makedirs(gene_table_dir, exist_ok=True)

    for flank in flank_sizes:

        start = gene_start - flank
        end = gene_end + flank

        final_start = max(0, start)
        final_end = min(chrom_len, end)

        actual_start_flank = gene_start - final_start
        actual_end_flank = final_end - gene_end

        cmd = [
            "klumpy",
            "alignment_plot",
            "--alignment_map", args.bam,
            "--reference", chrom,
            "--leftbound", str(final_start),
            "--rightbound", str(final_end),
            "--vertical_line_gaps",
            "--group_seqs",
            "--write_table", 
            "--annotation", args.annotation
        ]

        print("Running:", " ".join(cmd))

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print("ERROR running klumpy for gene:", gene, "flank:", flank)
            print(result.stderr)
            continue

        # move newest PDF
        pdf_files = glob.glob("*.pdf")
        newest_pdf = max(pdf_files, key=os.path.getctime)

        new_pdf = os.path.join(
            gene_plot_dir,
            f"{gene}_flank_{flank//1000}kb_alignment_plot.pdf"
        )
        shutil.move(newest_pdf, new_pdf)

        print("Saved plot:", new_pdf)

        # move newest table
        table_files = glob.glob("*.tsv")
        newest_table = max(table_files, key=os.path.getctime)

        new_table = os.path.join(
            gene_table_dir,
            f"{gene}_flank_{flank//1000}kb_table.tsv"
        )
        shutil.move(newest_table, new_table)

        print("Saved table:", new_table)

        summary.append({
            "focal_gene": gene,
            "requested_flank": flank,
            "actual_start_flank": actual_start_flank,
            "actual_end_flank": actual_end_flank
        })

summary_df = pd.DataFrame(summary)

summary_df.to_csv(
    args.summary,
    sep="\t",
    index=False
)

print("Finished generating alignment plots.")
