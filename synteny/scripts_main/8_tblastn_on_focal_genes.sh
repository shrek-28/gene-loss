#!/bin/bash

# Purpose:
# Run tblastn for all focal genes against all genomes
# and select best hit using:
#   1. highest query coverage
#   2. then highest bitscore

set -euo pipefail

# -------- ARGUMENT CHECK --------
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <FASTA_BASE> <GENOMES_DIR> <OUTPUT_DIR>"
    exit 1
fi

FASTA_BASE="$1"
GENOMES_DIR="$2"
OUTPUT_DIR="$3"

EVALUE=1e-10
MAX_TARGET_SEQS=20

mkdir -p "$OUTPUT_DIR"

# -------- LOOP OVER GENES --------
for GENE_DIR in "$FASTA_BASE"/*/; do
    GENE=$(basename "$GENE_DIR")
    echo "[INFO] Processing focal gene: $GENE"

    mkdir -p "$OUTPUT_DIR/$GENE"

    QUERY="$GENE_DIR/${GENE}.faa"
    if [ ! -f "$QUERY" ]; then
        echo "[WARNING] Missing FASTA: $QUERY"
        continue
    fi

    # -------- LOOP OVER GENOMES --------
    for DB in "$GENOMES_DIR"/*; do
        SPECIES=$(basename "$DB")

        OUTFILE="$OUTPUT_DIR/$GENE/${GENE}_vs_${SPECIES}.tsv"

        echo "[INFO] BLASTing $GENE vs $SPECIES"

        tblastn \
            -query "$QUERY" \
            -db "$DB/$SPECIES" \
            -evalue $EVALUE \
            -outfmt "6 qseqid sseqid sstart send evalue bitscore qcovs" \
            -max_target_seqs $MAX_TARGET_SEQS \
            -num_threads 8 \
            > "$OUTFILE"

    done
done

echo "[DONE] All focal gene tblastn searches completed."