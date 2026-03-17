#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pysam
import argparse
import os
import glob
import re

############################################
# ARGUMENTS
############################################
parser = argparse.ArgumentParser()
parser.add_argument("--table_root", required=True,
                    help="Root directory containing gene subdirectories with --write_table outputs")
parser.add_argument("--scan_alignments", required=True, help="scan_alignments TSV file")
parser.add_argument("--bam", required=True, help="BAM file (sorted and indexed)")
parser.add_argument("--region_table", required=True, help="TSV with focal_gene, chromosome, leftbound, rightbound")
parser.add_argument("--outdir", required=True, help="Output root directory for metrics")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

############################################
# LOAD DATA
############################################
scan_df = pd.read_csv(args.scan_alignments, sep="\t", comment="#")
region_df = pd.read_csv(args.region_table, sep="\t")
bam = pysam.AlignmentFile(args.bam, "rb")

############################################
# DISCOVER GENES
############################################
gene_dirs = [d for d in os.listdir(args.table_root)
             if os.path.isdir(os.path.join(args.table_root, d))]

############################################
# PROCESS EACH GENE
############################################
for gene in gene_dirs:
    gene_dir = os.path.join(args.table_root, gene)
    out_gene_dir = os.path.join(args.outdir, gene)
    os.makedirs(out_gene_dir, exist_ok=True)

    if gene not in region_df["focal_gene"].values:
        continue

    region = region_df[region_df["focal_gene"] == gene].iloc[0]
    chrom = region["chromosome"]
    start = int(region["leftbound"])
    end = int(region["rightbound"])
    region_len = end - start

    #########################################
    # COMMON BAM METRICS
    #########################################
    mapq_vals = []
    indel_count = 0
    for read in bam.fetch(chrom, start, end):
        mapq_vals.append(read.mapping_quality)
        if read.cigartuples:
            for op, length in read.cigartuples:
                if op in [1, 2]:  # Insertion/Deletion
                    indel_count += 1
    mapq_avg = np.mean(mapq_vals) if mapq_vals else 0

    coverage_vals = [p.n for p in bam.pileup(chrom, start, end)]
    coverage_mean = np.mean(coverage_vals) if coverage_vals else 0
    coverage_median = np.median(coverage_vals) if coverage_vals else 0
    coverage_sd = np.std(coverage_vals) if coverage_vals else 0

    common_metrics = {
        "chromosome": chrom,
        "gene_start": start,
        "gene_end": end,
        "region_length": region_len,
        "coverage_mean": coverage_mean,
        "coverage_median": coverage_median,
        "coverage_sd": coverage_sd,
        "mapq_avg": mapq_avg,
        "indel_count": indel_count
    }

    common_df = pd.DataFrame(common_metrics.items(), columns=["Metric", "Value"])
    common_df.to_csv(os.path.join(out_gene_dir, f"{gene}_common_metrics.tsv"),
                     sep="\t", index=False)

    #########################################
    # THRESHOLD METRICS
    #########################################
    table_files = glob.glob(os.path.join(gene_dir, "*_table.tsv"))
    threshold_rows = []

    for table in table_files:
        name = os.path.basename(table)
        m = re.search(r"flank_(\d+kb)", name)
        if not m:
            continue
        threshold = m.group(1)
        align_df = pd.read_csv(table, sep="\t")

        num_reads = len(align_df)
        num_groups = align_df["Group_Num"].nunique()
        group_sizes = align_df["Group_Num"].value_counts()
        largest_group_size = group_sizes.max()
        largest_group_fraction = largest_group_size / num_reads if num_reads else 0
        percent_aligned_avg = align_df["Percent_Aligned"].mean()
        clipped_start = align_df["Clipped_Start"].notna().sum()
        clipped_end = align_df["Clipped_End"].notna().sum()
        clipped_reads = clipped_start + clipped_end
        clipping_fraction = clipped_reads / num_reads if num_reads else 0
        alignment_span = (align_df["Position"] + align_df["Sequence_Length"]).max() - align_df["Position"].min()

        # Scan alignment contribution
        windows = scan_df[(scan_df["Reference_Seq"] == chrom) &
                          (scan_df["Start"] <= end) &
                          (scan_df["End"] >= start)]
        windows["Number_of_Groups"] = pd.to_numeric(windows["Number_of_Groups"], errors="coerce")

        if len(windows) == 0:
            S_scan = 1.0
        elif len(windows) <= 2:
            S_scan = 0.75
        elif len(windows) <= 5:
            S_scan = 0.4
        else:
            S_scan = 0.1

        # Normalized metrics
        EXPECTED_READS = 30
        S_reads = min(np.log1p(num_reads)/np.log1p(EXPECTED_READS),1)
        S_group = largest_group_fraction
        S_clip = max(1 - clipping_fraction, 0)
        S_align = min(percent_aligned_avg/100, 1)

        eps = 1e-6
        locus_score = ((S_reads+eps)*(S_group+eps)*(S_clip+eps)*(S_align+eps)*(S_scan+eps))**(1/5)

        threshold_rows.append({
            "threshold": threshold,
            "num_reads": num_reads,
            "num_groups": num_groups,
            "largest_group_fraction": largest_group_fraction,
            "clipping_fraction": clipping_fraction,
            "percent_aligned_avg": percent_aligned_avg,
            "alignment_span": alignment_span,
            "S_reads": S_reads,
            "S_group": S_group,
            "S_clip": S_clip,
            "S_align": S_align,
            "S_scan": S_scan,
            "locus_score": locus_score
        })

    thresh_df = pd.DataFrame(threshold_rows)
    thresh_df.to_csv(os.path.join(out_gene_dir, f"{gene}_threshold_metrics.tsv"),
                     sep="\t", index=False)

bam.close()