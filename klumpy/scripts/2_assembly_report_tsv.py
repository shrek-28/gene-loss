#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="Extract RefSeq chromosome mapping from NCBI assembly report"
)

parser.add_argument(
    "-i",
    "--assembly_report",
    required=True,
    help="Path to NCBI assembly_report.txt"
)

parser.add_argument(
    "-o",
    "--output",
    required=True,
    help="Output TSV mapping file"
)

args = parser.parse_args()

assembly_report = args.assembly_report
output_file = args.output

rows = []

with open(assembly_report) as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")

        sequence_name = fields[0]
        sequence_role = fields[1]
        refseq_accn = fields[6]

        # keep only chromosomes
        if sequence_role != "assembled-molecule":
            continue

        # remove version suffix
        refseq_clean = refseq_accn

        rows.append({
            "RefSeq": refseq_clean,
            "Chromosome": sequence_name
        })

df = pd.DataFrame(rows)

df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print(f"[DONE] Mapping written to {output_file}")