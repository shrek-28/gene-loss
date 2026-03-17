#!/usr/bin/env python3
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import pdist, squareform
from Bio import Phylo
import numpy as np
import argparse
from collections import defaultdict
import csv

# -------------------- Arguments --------------------
parser = argparse.ArgumentParser()
parser.add_argument("input_tree", help="Newick tree file")
parser.add_argument("output_file", help="Output figure file (PNG/PDF)")
parser.add_argument("--clusters", type=int, default=10)
parser.add_argument("--rep_file", help="Output CSV of representative species", default=None)
args = parser.parse_args()

# -------------------- Load tree --------------------
tree = Phylo.read(args.input_tree, "newick")
leaves = [leaf.name for leaf in tree.get_terminals()]
n = len(leaves)

# -------------------- Compute distance matrix --------------------
dm = np.zeros((n, n))
for i, l1 in enumerate(leaves):
    for j, l2 in enumerate(leaves):
        dm[i, j] = tree.distance(l1, l2)
dist_array = squareform(dm)

# -------------------- Hierarchical clustering --------------------
Z = linkage(dist_array, method='average')
clusters = fcluster(Z, t=args.clusters, criterion='maxclust')

# -------------------- Determine cut height --------------------
cut_height = sorted(Z[:,2], reverse=True)[args.clusters-1]

# -------------------- Plot dendrogram --------------------
plt.figure(figsize=(16,10))
dendrogram(
    Z,
    labels=leaves,
    leaf_rotation=90,
    leaf_font_size=10,
    color_threshold=cut_height  # auto-colors the 10 clusters
)
plt.axhline(y=cut_height, color='red', linestyle='--', linewidth=2, label=f'{args.clusters} clusters')
plt.title(f'Dendrogram with {args.clusters} clusters', fontsize=16)
plt.ylabel('Distance')
plt.xlabel('Leaf')
plt.legend()
plt.tight_layout()
plt.savefig(args.output_file, dpi=300)
print(f"Dendrogram saved to {args.output_file}")

# -------------------- Compute representative species --------------------
cluster_leaves = defaultdict(list)
for idx, cluster_id in enumerate(clusters):
    cluster_leaves[cluster_id].append(idx)

representatives = {}
for cluster_id, indices in cluster_leaves.items():
    if len(indices) == 1:
        representatives[cluster_id] = leaves[indices[0]]
    else:
        intra_dist = dm[np.ix_(indices, indices)]
        sum_dist = intra_dist.sum(axis=1)
        rep_idx = indices[np.argmin(sum_dist)]
        representatives[cluster_id] = leaves[rep_idx]

# Print representative species
print("\nRepresentative species per cluster:")
for cluster_id in sorted(representatives.keys()):
    print(f"Cluster {cluster_id}: {representatives[cluster_id]}")

# Save to CSV if requested
if args.rep_file:
    with open(args.rep_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Cluster', 'Representative'])
        for cluster_id in sorted(representatives.keys()):
            writer.writerow([cluster_id, representatives[cluster_id]])
    print(f"Representative species saved to {args.rep_file}")