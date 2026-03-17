## STEP 1
## For Step 1, download all genomes from HGNC. 
## Link: https://www.genenames.org/download/ 
## download the TXT file by 'save page as' from complete set new TSV
## saved file in .txt format.

## conversion of txt to tsv format 
## replace filenames with what you have saved as
mv filename.txt filename.tsv

## generating a file with col headers and col numbers
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) print i, $i}' file.tsv > filename.txt

## col 19 -> entrez ID/gene ID
## extracting col 19 values: 
awk -F'\t' 'NR>1 {print $19}' geneids.tsv | sort -u > geneids_list.txt

## geneids_list.txt is the list of Entrez IDs of all genes from HGNC