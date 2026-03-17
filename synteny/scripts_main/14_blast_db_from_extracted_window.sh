#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# -----------------------------
# ARGUMENT PARSER
# -----------------------------
usage() {
    echo "Usage: $0 --extracted_root <path> --db_root <path>"
    exit 1
}

EXTRACTED_ROOT=""
DB_ROOT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --extracted_root)
            EXTRACTED_ROOT="$2"
            shift 2
            ;;
        --db_root)
            DB_ROOT="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            usage
            ;;
    esac
done

if [[ -z "$EXTRACTED_ROOT" || -z "$DB_ROOT" ]]; then
    usage
fi

mkdir -p "$DB_ROOT"

echo "[INFO] Creating BLAST databases from extracted regions"

# -----------------------------
# LOOP OVER FOCAL GENES
# -----------------------------
for gene_dir in "$EXTRACTED_ROOT"/*; do
    [[ ! -d "$gene_dir" ]] && continue

    focal_gene=$(basename "$gene_dir")
    echo "[INFO] Processing focal gene: $focal_gene"

    gene_db_dir="$DB_ROOT/$focal_gene"
    mkdir -p "$gene_db_dir"

    # -----------------------------
    # LOOP OVER SPECIES
    # -----------------------------
    for species_dir in "$gene_dir"/*; do
        [[ ! -d "$species_dir" ]] && continue

        species=$(basename "$species_dir")

        fasta_files=( "$species_dir"/*.fna )

        if [[ ${#fasta_files[@]} -eq 0 ]]; then
            echo "[WARNING] No FASTA found for $focal_gene | $species"
            continue
        fi

        out_species_db_dir="$gene_db_dir/$species"
        mkdir -p "$out_species_db_dir"

        db_prefix="$out_species_db_dir/${focal_gene}_${species}_region_db"

        echo "  [RUNNING] makeblastdb for $focal_gene | $species"

        makeblastdb \
            -in "${fasta_files[0]}" \
            -dbtype nucl \
            -parse_seqids \
            -out "$db_prefix" \
            > /dev/null

        echo "  [OK] DB created"
    done
done

echo "[DONE] All BLAST databases created."