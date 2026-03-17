#!/usr/bin/env python3
"""
Clean FASTA headers to keep only the XP protein_id.
Usage:
    python clean_fasta_headers.py input.faa output_clean.faa
"""

import sys
from Bio import SeqIO
import re

if len(sys.argv) != 3:
    print("Usage: python clean_fasta_headers.py input.faa output_clean.faa")
    sys.exit(1)

input_faa = sys.argv[1]
output_faa = sys.argv[2]

with open(output_faa, "w") as out_handle:
    for record in SeqIO.parse(input_faa, "fasta"):
        # Extract XP protein ID from header
        m = re.search(r'^(XP_\d+\.\d+)', record.description)
        if m:
            xp_id = m.group(1)
        else:
            continue  

        record.id = xp_id
        record.description = ""  # remove everything else
        SeqIO.write(record, out_handle, "fasta")

print(f"[OK] Cleaned headers written to {output_faa}")