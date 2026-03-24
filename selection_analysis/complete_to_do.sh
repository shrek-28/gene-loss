## STEP 0: DOWNLOAD GTF AND GENOME FASTA FILES OF ALL DESIRED SPECIES

## STEP 1: OBTAINING CDS 
mkdir -p data/cds_sequences
while read species; do
    ./scripts/cds_extraction.sh \
        data/gtf/${species}.gff \
        data/fasta/${species}.fna \
        data/cds_sequences/${species}_cds.fa
done < species_list.txt