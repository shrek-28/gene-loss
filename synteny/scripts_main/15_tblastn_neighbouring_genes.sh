#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# -----------------------------
# ARGUMENT PARSER
# -----------------------------
usage() {
    echo "Usage: $0 --fa_dir <path> --db_dir <path> --out_dir <path> [--threads <int>]"
    exit 1
}

FA_DIR=""
DB_DIR=""
OUT_DIR=""
THREADS=4

while [[ $# -gt 0 ]]; do
    case $1 in
        --fa_dir)
            FA_DIR="$2"
            shift 2
            ;;
        --db_dir)
            DB_DIR="$2"
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

if [[ -z "$FA_DIR" || -z "$DB_DIR" || -z "$OUT_DIR" ]]; then
    usage
fi

mkdir -p "$OUT_DIR"

echo "[INFO] Starting tblastn for neighbouring genes"

# -----------------------------
# LOOP OVER FOCAL GENES
# -----------------------------
for focal_gene_dir in "$DB_DIR"/*; do
    [[ ! -d "$focal_gene_dir" ]] && continue

    focal_gene=$(basename "$focal_gene_dir")
    echo "[INFO] Processing focal gene: $focal_gene"

    fasta_dir="$FA_DIR/$focal_gene"
    [[ ! -d "$fasta_dir" ]] && {
        echo "[WARNING] Missing FASTA dir for $focal_gene"
        continue
    }

    # -----------------------------
    # LOOP OVER NEIGHBOUR GENES
    # -----------------------------
    for faa_file in "$fasta_dir"/*.faa; do

        file_base=$(basename "$faa_file")

        # Skip focal gene itself
        if [[ "$file_base" == "$focal_gene.faa" ]]; then
            continue
        fi

        neigh_gene=$(basename "$faa_file" .faa)
        echo "  [RUNNING] neighbour: $neigh_gene"

        neigh_out_dir="$OUT_DIR/$focal_gene/$neigh_gene"
        mkdir -p "$neigh_out_dir"

        # -----------------------------
        # LOOP OVER SPECIES DBs
        # -----------------------------
        for species_dir in "$focal_gene_dir"/*; do

            [[ ! -d "$species_dir" ]] && continue

            species=$(basename "$species_dir")
            db_prefix="$species_dir/${focal_gene}_${species}_region_db"

            if [[ ! -f "${db_prefix}.nin" ]]; then
                echo "    [SKIP] Missing DB for $species"
                continue
            fi

            out_file="$neigh_out_dir/${focal_gene}_${neigh_gene}_vs_${species}.tblastn.tsv"

            tblastn \
                -query "$faa_file" \
                -db "$db_prefix" \
                -out "$out_file" \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
                -num_threads "$THREADS"

            echo "    [DONE] $neigh_gene vs $species"
        done
    done
done

echo "[DONE] All neighbouring gene tblastn completed."