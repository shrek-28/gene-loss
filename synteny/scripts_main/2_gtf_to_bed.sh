#!/bin/bash
# Usage: ./gtf_to_bed.sh input.gtf output.bed

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.gtf output.bed"
    exit 1
fi

INPUT="$1"
OUTPUT="$2"

# match($0, /GeneID:([0-9]+)/, arr) # was used when Entrez ID is extracted 
## change as per need

awk 'BEGIN{OFS="\t"; print "chrom\tstart\tend\tgeneID\tnone\tstrand"}
!/^#/ {
    if($3=="gene"){
        chrom=$1
        start=$4-1       # 0-based BED
        end=$5
        strand=$7
        match($0, /gene "([^"]+)"/, arr)
        gene_id=arr[1]
        print chrom, start, end, gene_id, ".", strand
    }
}' "$INPUT" > "$OUTPUT"
