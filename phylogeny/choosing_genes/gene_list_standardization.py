# Input file path
input_file = "phylogeny/choosing_genes/gene_list.txt"

# Read, capitalize, and store in list
with open(input_file, "r") as f:
    gene_list = [line.strip().upper() for line in f if line.strip()]

# Remove duplicates
unique_genes = list(set(gene_list))

# Sort alphabetically
unique_genes.sort()

# Print result
print(unique_genes)

# Optional: write to a new file
with open("phylogeny/choosing_genes/unique_sorted_genes.txt", "w") as out:
    for gene in unique_genes:
        out.write(gene + "\n")