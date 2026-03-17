#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 --neigh_dir <dir> --extracted_dir <dir> --out_file <file>"
    exit 1
}

NEIGH_DIR=""
EXTRACTED_DIR=""
OUT_FILE=""

# -----------------------------
# ARGUMENT PARSER
# -----------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --neigh_dir)
            NEIGH_DIR="$2"
            shift 2
            ;;
        --extracted_dir)
            EXTRACTED_DIR="$2"
            shift 2
            ;;
        --out_file)
            OUT_FILE="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            usage
            ;;
    esac
done

if [[ -z "$NEIGH_DIR" || -z "$EXTRACTED_DIR" || -z "$OUT_FILE" ]]; then
    usage
fi

# -----------------------------
# HEADER
# -----------------------------
echo -e "focal_gene\tgene_id\tgene_type\tseq_len\tblast_task" > "$OUT_FILE"

# -----------------------------
# LOOP OVER NEIGHBOUR FILES
# -----------------------------
for neigh_file in "$NEIGH_DIR"/*_neighbour_summary.tsv; do

    focal_gene=$(basename "$neigh_file" _neighbour_summary.tsv)

    tail -n +2 "$neigh_file" | while IFS=$'\t' read -r gene_id gene_type _; do

        # Skip protein coding
        if [[ "$gene_type" == "protein_coding" ]]; then
            continue
        fi

        query_fa="$EXTRACTED_DIR/$focal_gene/$gene_id.fa"

        if [[ ! -f "$query_fa" ]]; then
            continue
        fi

        # Sequence length (ignore header)
        seq_len=$(grep -v ">" "$query_fa" | tr -d '\n' | wc -c)

        # Choose BLAST task
        if [[ "$gene_type" == "snoRNA" || "$seq_len" -lt 200 ]]; then
            blast_task="blastn-short"
        else
            blast_task="dc-megablast"
        fi

        echo -e "${focal_gene}\t${gene_id}\t${gene_type}\t${seq_len}\t${blast_task}" >> "$OUT_FILE"

    done
done

echo "[DONE] Summary TSV created at $OUT_FILE"