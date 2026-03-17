#!/usr/bin/env python3
import os
import re
import argparse
from collections import defaultdict
from Bio import SeqIO

# --------------------------
# Argument parser
# --------------------------
parser = argparse.ArgumentParser(
    description="Extract N upstream and downstream gene neighbors with gene type and longest protein."
)

parser.add_argument("--bed", required=True, help="Path to genes.bed file")
parser.add_argument("--gff", required=True, help="Path to GTF/GFF file")
parser.add_argument("--proteins", required=True, help="Path to protein FASTA (.faa)")
parser.add_argument("--gene_ids", required=True, help="Path to focal gene ID list")
parser.add_argument("--output_dir", required=True, help="Directory to write output TSVs")
parser.add_argument("--neighbors", type=int, default=5,
                    help="Number of upstream and downstream neighbors (default=5)")

args = parser.parse_args()

BED_FILE = args.bed
GFF_FILE = args.gff
PROTEIN_FASTA = args.proteins
GENE_IDS = args.gene_ids
OUTPUT_DIR = args.output_dir
N = args.neighbors

os.makedirs(OUTPUT_DIR, exist_ok=True)

# --------------------------
# Step 1: Read BED file
# --------------------------
genes_by_chr = defaultdict(list)
gene_info = dict()

with open(BED_FILE) as f:
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) < 6:
            continue
        chrom, start, end, gene, score, strand = fields[:6]
        genes_by_chr[chrom].append(gene)
        gene_info[gene] = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand
        }

# --------------------------
# Step 2: Parse GTF/GFF
# --------------------------
gene_type_map = dict()
gene_proteins = defaultdict(list)

with open(GFF_FILE) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue

        feature_type = fields[2]
        attrs = fields[8]

        # Extract Entrez GeneID using regex
        match = re.search(r'gene "([^"]+)"', attrs)        
        gene_match = match.group(1) if match else None

        if feature_type == "gene" and gene_match:
            type_match = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_biotype") or attr.startswith("gene_type"):
                    type_match = attr.split('"')[1]

            gene_type_map[gene_match] = type_match if type_match else "unknown"

        elif feature_type == "CDS" and gene_match:
            protein_match = None
            for attr in attrs.split(";"):
                attr = attr.strip()
                if attr.startswith("protein_id"):
                    protein_match = attr.split('"')[1]

            if protein_match:
                gene_proteins[gene_match].append(protein_match)

# --------------------------
# Step 3: Protein lengths
# --------------------------
prot_len = dict()
for rec in SeqIO.parse(PROTEIN_FASTA, "fasta"):
    prot_len[rec.id] = len(rec.seq)

# --------------------------
# Step 4: Longest protein per gene
# --------------------------
longest_protein = dict()

for gene, prots in gene_proteins.items():
    if prots:
        best = max(prots, key=lambda x: prot_len.get(x, 0))
        longest_protein[gene] = best

# --------------------------
# Step 5: Read focal genes
# --------------------------
with open(GENE_IDS) as f:
    focal_genes = [line.strip() for line in f if line.strip()]

# --------------------------
# Step 6: Extract neighbors
# --------------------------
for focal in focal_genes:

    chrom = gene_info.get(focal, {}).get("chrom")
    if not chrom or chrom not in genes_by_chr:
        print(f"[WARNING] Focal gene {focal} not found in BED. Skipping.")
        continue

    genes_on_chr = genes_by_chr[chrom]
    idx = genes_on_chr.index(focal)

    start = max(0, idx - N)
    end = min(len(genes_on_chr), idx + N + 1)
    neighborhood = genes_on_chr[start:end]

    out_file = os.path.join(OUTPUT_DIR, f"{focal}_neighbour_summary.tsv")

    with open(out_file, "w") as out:
        out.write("gene_id\tgene_type\tlongest_protein\n")
        for g in neighborhood:
            gtype = gene_type_map.get(g, "unknown")
            prot = longest_protein.get(g, "NA") if gtype == "protein_coding" else gtype
            out.write(f"{g}\t{gtype}\t{prot}\n")

print(f"Neighbour summary TSVs written to {OUTPUT_DIR}")
