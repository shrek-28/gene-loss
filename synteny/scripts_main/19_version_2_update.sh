#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

usage() {
    echo "Usage: $0 --noncoding_dir <dir> --db_dir <dir> --tsv <file> --out_dir <dir> [--threads <n>]"
    exit 1
}

NONCODING_DIR=""
DB_DIR=""
TSV=""
OUT_DIR=""
THREADS=4

# -----------------------------
# ARGUMENT PARSER
# -----------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --noncoding_dir)
            NONCODING_DIR="$2"
            shift 2
            ;;
        --db_dir)
            DB_DIR="$2"
            shift 2
            ;;
        --tsv)
            TSV="$2"
            shift 2
            ;;
        --out_dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            usage
            ;;
    esac
done

if [[ -z "$NONCODING_DIR" || -z "$DB_DIR" || -z "$TSV" || -z "$OUT_DIR" ]]; then
    usage
fi

mkdir -p "$OUT_DIR"

# -----------------------------
# LOOP THROUGH TSV
# -----------------------------
tail -n +2 "$TSV" | while IFS=$'\t' read -r focal_gene gene_id gene_type seq_len blast_task; do

    query_fa="$NONCODING_DIR/$focal_gene/$gene_id.fa"
    [[ ! -f "$query_fa" ]] && { echo "[WARN] Missing query: $query_fa"; continue; }

    gene_db_dir="$DB_DIR/$focal_gene"
    [[ ! -d "$gene_db_dir" ]] && { echo "[WARN] Missing DB dir: $gene_db_dir"; continue; }

    echo "[INFO] BLASTing $gene_id ($blast_task)"

    # -----------------------------
    # LOOP OVER SPECIES DBs
    # -----------------------------
    for species_dir in "$gene_db_dir"/*; do

        species=$(basename "$species_dir")
        db_prefix="$species_dir/${focal_gene}_${species}_region_db"

        [[ ! -f "${db_prefix}.nin" ]] && continue

        outdir="$OUT_DIR/$focal_gene/$gene_id/$species"
        mkdir -p "$outdir"

        outfile="$outdir/${gene_id}_vs_${species}.blastn.tsv"

        blastn \
            -task "$blast_task" \
            -query "$query_fa" \
            -db "$db_prefix" \
            -out "$outfile" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -evalue 1e-5 \
            -max_target_seqs 10 \
            -num_threads "$THREADS"

        if [[ ! -s "$outfile" ]]; then
            echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs" > "$outfile"
        fi

    done
done

echo "[DONE] All BLAST searches completed"