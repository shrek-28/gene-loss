## ALL SCRIPTS SAVED IN SUBFOLDER "scripts" AND ALL DATA PRESENT IN SUBFOLDER "data"

## Script 1 - making gene list 
## script to be executed line by line 
## file saved as data/unique_sorted_genes.txt

## step 1a - downloading human genome (reference genome)
# genome FASTA, GTF and protein FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz > data/human_genome_files/GCF_000001405.40_GRCh38.p14_genomic.fna
gunzip GCF_000001405.40_GRCh38.p14_genomic.gtf.gz > data/human_genome_files/GCF_000001405.40_GRCh38.p14_genomic.gtf
gunzip GCF_000001405.40_GRCh38.p14_protein.faa.gz > data/human_genome_files/GCF_000001405.40_GRCh38.p14_protein.faa

## step 2 - GTF to BED conversion of reference genome 
./scripts/2_gtf_to_bed.sh data/human_genome_files/GCF_000001405.40_GRCh38.p14_genomic.gtf data/human_genome_files/human_bed_with_gene_symbol.bed
## NOTE: CHANGE SCRIPT IF ANY OTHER IDENTIFIER OTHER THAN GENE SYMBOL IS USED

## step 3 - cleaning of FASTA headers
python3 scripts/3_cleaning_fasta_headers.py data/human_genome_files/GCF_000001405.40_GRCh38.p14_protein.faa data/human_clean_protein_fasta.faa
## step 3a - only to retain XP
python3 scripts/3a_cleaning_fasta_headers_only_xp.py data/human_genome_files/GCF_000001405.40_GRCh38.p14_protein.faa data/clean_human_only_XP.faa

## step 4 - extraction of neighbours
python3 scripts/4_extract_neighbours.py data/human_genome_files/human_bed_with_gene_symbol.bed data/unique_sorted_genes.txt data/neighbors_output

## step 5 - extract neighbourhood coordinates 
python3 scripts/5_extract_neighbourhood_coords.py data/neighbours_output data/human_genome_files/human_bed_with_gene_symbol.bed data/neighbour_summary 

## step 6 - generating neighbour TSV
## NOTE: CHANGE SCRIPT IF ANY OTHER IDENTIFIER OTHER THAN GENE SYMBOL IS USED
python3 scripts/6_generating_neighbour_tsv.py --bed data/neighbours_output data/human_genome_files/human_bed_with_gene_symbol.bed --gtf data/human_genome_files/GCF_000001405.40_GRCh38.p14_genomic.gtf --proteins data/human_genome_files/GCF_000001405.40_GRCh38.p14_protein.faa --gene_ids data/unique_sorted_genes.txt --output_dir data/neighbour_summary

## step 7 - obtain protein sequence 
python3 scripts/7_obtaining_protein_sequences.py --proteins data/human_clean_protein_fasta.faa --output_dir data/results_per_gene_fasta_py

## step 7a - generation of blast database 
## list of all species saved as species_list.txt
## prereq: download all genome files into one folder - blast_database
## run the below code
for f in blast_database/genome_FASTA/*.fna; do
    base=$(basename "$f" .fna)

    mkdir -p blastdb/"$base"

    makeblastdb \
        -in "$f" \
        -dbtype nucl \
        -parse_seqids \
        -out blastdb/"$base"/"$base"
done

## step 8 - run tblastn on focal genes
./scripts/8_tblastn_on_focal_genes.sh data/results_per_gene_fasta_py data/blastdb data/results_per_gene_tblastn

## step 9 - obtain best hits 
python3 scripts/9_obtain_best_hits.py -t data/results_per_gene_tblastn -n data/neighbour_summary -o data/focal_gene_best_hits

## step 10 - obtaining best hits with status 
python3 scripts/10_obtain_best_hits_with_status.py -i data/focal_gene_best_hits -o data/focal_gene_best_hits_with_status

## step 11 - extracting window length
python3 scripts/11_get_extract_window_length.py -b data/focal_gene_best_hits -n data/neighbour_summary -o data/extraction_windows

## step 12 - validation of windows
## do only for reference genome - human in this case 
python3 scripts/12_validate_windows.py -n data/neighbour_summary -w data/extraction_windows -s Homo_sapiens

## step 13 - extract window from genome 
python3 scripts/13_extract_window_from_genome.py --window_root data/extracted_windows --genome_dir data/genomes_dir --output_root data/extracted_regions

## step 14 - making blastdb from extracted region 
./scripts/14_make_database_from_extracted_window.sh --extracted_root data/extracted_regions --db_root data/extracted_regions_db

## script 15 - tblastn neighbouring genes search 
./scripts/15_tblastn_neighbouring_genes.sh --fa_dir data/results_per_gene_fasta_py --db_dir data/extracted_regions_db --out_dir data/tblastn_results_neighbours --threads 8

## script 16 - obtain neighbour best hit 
python3 scripts/16_obtain_neighbour_best_hit.py --root_dir data/tblastn_results_neighbours

## script 17 - obtain fasta of non protein coding genes
./scripts/17_extract_nonproteingcoding_fasta.sh --genome data/genomes_dir/Homo_sapiens.fna --neigh_dir data/neighbour_summary --genes_bed data/human_genome_files/human_bed_with_gene_symbol.bed --outdir data/non-proteincodinggenes

## script 18 - checking which blast 
./scripts/18_check_which_blast.sh --neigh_dir data/neighbour_summary --extracted_dir data/non-proteincodinggenes --out_file data/check_blast_nonproteincoding.tsv

## script 19 - blast non protein coding 
./scripts/19_version2update.sh --noncoding_dir data/non-proteincodinggenes --db_dir data/extracted_regions_db --tsv data/check_blast_nonproteincoding.tsv --out_dir data/blast_non_proteincoding --threads 8

## script 20 - best hits non protein coding 
python3 scripts/20_obtain_besthits_for_nonproteincoding.py --root_dir data/blast_non_proteincoding

## script 21 - add genomic coordinates for protein coding 
python3 scripts/21_add_genomic_coordinates_protein_coding.py -i data/tblastn_results_neighbours -o data/neighbouring_genes_genomic_coordinates

## script 22 - add genomic coordinates non protein coding
python3 scripts/22_add_genomic_coordinates_nonprotein_coding.py -i data/blast_non_proteincoding -o data/neighbouring_genes_genomic_coordinates

## script 23 - adding focal genes 
python3 scripts/23_add_focal_genes.py -i data/focal_gene_best_hits_with_status -o data/neighbouring_genes_genomic_coordinates

## script 24 - generation of synteny tables 
python3 scripts/24_synteny_table_creation.py --coord_dir data/neighbour_summary --gene_coord_root data/neighbouring_genes_genomic_coordinates --species_list data/species_list.txt --out_root data/synteny_tables

## script 24a - klumpy synteny tables
python3 scripts/synteny_table_for_klumpy.py --coord_dir data/neighbour_summary --gene_coord_root data/neighbouring_genes_genomic_coordinates --species_list data/species_list.txt --out_root data/synteny_tables_for_klumpy

## script 25 - gene presence absence plots 
python3 scripts/25_gene_presence_absence_plots.py --neigh_coord_dir data/neighbour_summary --coord_root data/neighbouring_genes_genomic_coordinates --species_list data/species_list.txt --out_dir data/genepresencematrix

## script 26 - synteny metrics 
python3 scripts/26_synteny_metrics.py --table_dir data/synteny_tables --ref_dir data/neighbour_summary --out_dir data/metrics

## script 27 - generation of latex report 
python3 scripts/27_generate_latex_report.py --presence_dir data/genepresencematrix --metric_dir data/metrics --out_tex data/synteny_report.tex
pdflatex data/synteny_report.tex
pdflatex data/synteny_report.tex 
## run twice - imp.
