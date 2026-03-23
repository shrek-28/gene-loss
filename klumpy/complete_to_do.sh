## DOWNLOADING SEQUENCE FROM NCBI SRA 

# installation of NCBI SRA toolkit
mkdir -p ~/software
cd ~/software
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xvzf sratoolkit.current-ubuntu64.tar.gz

# temporary addition to PATH
export PATH=$HOME/software/sratoolkit.3.3.0-ubuntu64/bin:$PATH

# checks 
prefetch --version
fasterq-dump --version

## --------------------------------------------------------------------------------------------------------------------------------------------
## --------------------------------------------------------------------------------------------------------------------------------------------

## downloading the SRA file 

## getting SRR id
# bos taurus
esearch -db sra -query SRX24621186 | efetch -format runinfo | cut -d ',' -f 1
## RESULT: SRR29097113
# pseudorca crassidens
esearch -db sra -query SRX25433668 | efetch -format runinfo | cut -d ',' -f 1
## RESULT: SRR29940089

## installation
# bos taurus
prefetch --max-size 100G SRR29097113
## download is done with increased maximum download size 
# pseudorca crassidens
prefetch --max-size 100G SRR29940089
nohup prefetch SRR29940089 > prefetch.log 2>&1 &

## downloading using aspera
 /home/ceglab27/.aspera/connect/bin/ascp -i /home/ceglab27/.aspera/connect/etc/asperaweb_id_dsa.openssh \
     -QT -l 300m -P33001 \
    era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR299/089/SRR29940089/SRR29940089_subreads.fastq.gz \
     .
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

## fasterq dump and zipping 
# fasterq dump 
nohup fastq-dump SRR29097113 --split-spot --gzip --readids > fastqdump.log 2>&1 &

## supplementary 
vdb-dump --info SRR29097113 
## checks information to see how big the fastq file will be (expected 10-25 GB)

## post-fastq dump check of length 
zcat SRR29097113.fastq.gz | wc -l
## SEQ : 4,309,512 * 4 should be the value
17238048

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

## minimap2 installation 

#installation 
git clone https://github.com/lh3/minimap2
cd minimap2
make

# within minimap2 dir
# testing whether installation has happened
./minimap2 --help

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# alignment using minimap2 

## minimap2
nohup bash -c "minimap2/minimap2 -ax map-hifi -t 16 Bos_taurus.fna SRR29097113.fastq.gz | samtools sort -@ 16 -o aln.bam" > minimap2.log 2>&1 &
nohup bash -c "minimap2/minimap2 -ax map-pb -t 16 GCF_039906515.1_mPseCra1.hap1_genomic.fna SRR29940089_subreads.fastq.gz | samtools sort -@ 16 -o Pseudorca_crassidens.bam" > minimap2.log 2>&1 &

## checking if it is sorted
samtools view -H aln.bam | grep SO
## should give SO:coordinate

## indexing 
samtools index aln.bam

## checking chromosome patterns in bos taurus
samtools idxstats aln.bam

## CHECKING 
samtools view -H Pseudorca_crassidens.bam | less -S
## should contain @SQ lines 

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

## klumpy installation
pip install klumpy
## (samtools should be already present in path)

## bos taurus gtf download 
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_genomic.gtf.gz
gunzip -c GCF_002263795.3_ARS-UCD2.0_genomic.gtf.gz > Bos_taurus.gtf

# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# complete till step 23 of synteny protocol

## STEP 1 - REFINED SYNTENY TABLE CREATION 
python3 scripts/synteny_table_for_klumpy.py --coord_dir data/neighbour_summary --gene_coord_root data/neighbouring_genes_genomic_coordinates --species_list data/species_list_klumpy.txt --out_root data/synteny_tables_for_klumpy

## STEP1a - standardization of chromosome names
## bos taurus assembly report download 
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0/GCF_002263795.3_ARS-UCD2.0_assembly_report.txt
## pseudorca crassidens report download
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/039/906/515/GCF_039906515.1_mPseCra1.hap1/GCF_039906515.1_mPseCra1.hap1_assembly_report.txt
## parsing 
python3 scripts/assembly_report_tsv.py --assembly_report data/GCF_002263795.3_ARS-UCD2.0_assembly_report.txt --output data/bos_taurus_chr.tsv
python3 scripts/assembly_report_tsv.py --assembly_report data/GCF_039906515.1_mPseCra1.hap1_assembly_report.txt --output data/pseudorca_crassidens_chr.tsv
## replacement
python3 scripts/mapping_chr.py --input_dir data/synteny_tables_for_klumpy --output_dir data/mapped_synteny_tables --map Bos_taurus data/bos_taurus_chr.tsv --map Pseudorca_crassidens data/pseudorca_crassidens_chr.tsv

## STEP 2 - SPECIES WISE TABLES 
python3 scripts/species_wise_tables.py --input_dir data/mapped_synteny_tables --output_dir data/species_wise_klumpy_tables

## STEP3 - SCAN ALIGNMENTS = STEP INDEPENDENT OF PREVIOUS STEPS
python3 scripts/scan_alignments.py --species Bos_taurus --bam data/aln.bam --annotation data/Bos_taurus.gtf --threads 16 --outdir data/Bos_taurus_klumpy_results --flag_excess_groups 
## nohup
nohup python3 scripts/scan_alignments.py \
--species Bos_taurus \
--bam data/aln.bam \
--annotation data/Bos_taurus.gtf \
--threads 16 \
--outdir data/Bos_taurus_klumpy_results \
--flag_excess_groups \
> scan_alignments.log 2>&1 &
## pseudorca crassidens
nohup python3 scripts/scan_alignments.py \
--species Pseudorca_crassidens \
--bam data/Pseudorca_crassidens.bam \
--annotation data/Pseudorca_crassidens.gtf \
--threads 16 \
--outdir data/Pseudorca_crassidens_klumpy_results \
--flag_excess_groups \
> scan_alignments.log 2>&1 &

## STEP 4 - ALIGNMENT PLOTS 
python3 scripts/alignment_plots_progressive_new.py --bam data/Bos_taurus_files/aln.bam --tsv data/species_wise_klumpy_tables/Bos_taurus.tsv --outdir data/Bos_taurus_alignment_plots_progressive_new --table_dir data/Bos_taurus_alignment_tables_progressive_new --summary data/flanking_regions_table_progressive_new.tsv --annotation 
python3 scripts/alignment_plots.py --bam data/Pseudorca_crassidens.bam --tsv data/species_wise_klumpy_tables/Pseudorca_crassidens.tsv --outdir data/Pseudorca_crassidens_alignment_plots --table_dir data/Pseudorca_crassidens_alignment_tables --summary data/Pseudorca_crassidens_flanking_table.tsv --annotation data/Pseudorca_crassidens.gtf

## STEP 5 - ALIGNMENT METRICS
python3 scripts/alignment_metrics_progressive.py --table_root data/Bos_taurus_alignment_tables_progressive --scan_alignments data/scan_alignments_results/Bos_taurus_Candidate_Regions.tsv --bam data/Bos_taurus_files/aln.bam --region_table data/species_wise_klumpy_tables/Bos_taurus.tsv --outdir data/Bos_taurus_alignment_metrics_progressive

## STEP 6 - LATEX FILE CREATION
python3 scripts/latex_file_creation.py --plots_dir data/Bos_taurus_alignment_plots_progressive_new --metrics_dir data/Bos_taurus_alignment_metrics_progressive --flanking_tsv data/flanking_regions_table_progressive.tsv --output_tex Bos_taurus_progressive.tex --organism_name "Bos taurus" --output_pdf Bos_taurus_final_klumpy.pdf