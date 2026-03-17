#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 --genome <genome_fasta> --neigh_dir <dir> --genes_bed <file> --outdir <dir>"
    exit 1
}

GENOME=""
NEIGH_DIR=""
GENES_BED=""
OUTDIR=""

# -----------------------------
# ARGUMENT PARSER
# -----------------------------
while [[ $# -gt 0 ]]; do
    case $1 in
        --genome)
            GENOME="$2"
            shift 2
            ;;
        --neigh_dir)
            NEIGH_DIR="$2"
            shift 2
            ;;
        --genes_bed)
            GENES_BED="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1"
            usage
            ;;
    esac
done

if [[ -z "$GENOME" || -z "$NEIGH_DIR" || -z "$GENES_BED" || -z "$OUTDIR" ]]; then
    usage
fi

# -----------------------------
# DEPENDENCY CHECK
# -----------------------------
if ! command -v samtools &> /dev/null; then
    echo "[ERROR] samtools not found. Please install it first."
    exit 1
fi

# -----------------------------
# INDEX GENOME
# -----------------------------
if [[ ! -f "${GENOME}.fai" ]]; then
    echo "[INFO] Indexing genome with samtools faidx"
    samtools faidx "$GENOME"
fi

mkdir -p "$OUTDIR"

# -----------------------------
# PROCESS NEIGHBOUR FILES
# -----------------------------
for neigh_file in "${NEIGH_DIR}"/*_neighbour_summary.tsv; do

    focal_gene=$(basename "$neigh_file" _neighbour_summary.tsv)
    echo "[INFO] Processing focal gene: $focal_gene"

    focal_out="${OUTDIR}/${focal_gene}"
    mkdir -p "$focal_out"

    tail -n +2 "$neigh_file" | while IFS=$'\t' read -r gene_id gene_type rest; do

        # STRICT filter
        if [[ "$gene_type" == "protein_coding" ]]; then
            continue
        fi

        bed_line=$(awk -v g="$gene_id" '$4==g {print; exit}' "$GENES_BED")
        if [[ -z "$bed_line" ]]; then
            echo "  [WARN] Coordinates not found for $gene_id"
            continue
        fi

        chrom=$(echo "$bed_line" | cut -f1)
        bed_start=$(echo "$bed_line" | cut -f2)
        bed_end=$(echo "$bed_line" | cut -f3)

        # BED (0-based) → samtools (1-based)
        faidx_start=$((bed_start + 1))
        faidx_end=$bed_end

        out_fa="${focal_out}/${gene_id}.fa"

        echo "  [EXTRACT] $gene_id ($chrom:${faidx_start}-${faidx_end})"

        samtools faidx "$GENOME" "${chrom}:${faidx_start}-${faidx_end}" > "$out_fa"

    done
done

echo "[DONE] Only non-protein-coding neighbouring genes extracted"