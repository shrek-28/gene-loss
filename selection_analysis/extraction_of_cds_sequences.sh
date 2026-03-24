#!/bin/bash

# Usage:
# bash extract_cds_gff3.sh input.gff3 genome.fasta output_cds.fa

gff_file=$1
genome_fasta=$2
output_cds=$3

echo "Processing $gff_file..."

# Step 1: keeping only CDS features
awk '$3=="CDS"' "$gff_file" > tmp_cds.gff3

# Step 2: Extract CDS 
gffread tmp_cds.gff3 -g "$genome_fasta" -x "$output_cds"

# Cleanup
rm tmp_cds.gff3

echo "[DONE] CDS written to $output_cds"