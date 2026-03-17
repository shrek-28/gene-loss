#!/usr/bin/env python3
import os
import glob
import pandas as pd
import argparse

# --------------------------------------------------
# Argument parser
# --------------------------------------------------
parser = argparse.ArgumentParser(
    description="Generate LaTeX synteny verification report."
)

parser.add_argument(
    "--base_dir",
    required=False,
    help="Base directory of the synteny pipeline data"
)

parser.add_argument(
    "--presence_dir",
    required=False,
    help="Directory containing *_gene_presence.png files"
)

parser.add_argument(
    "--metric_dir",
    required=False,
    help="Directory containing *_metrics.tsv files"
)

parser.add_argument(
    "--out_tex",
    required=True,
    help="Output LaTeX file"
)

args = parser.parse_args()

# --------------------------------------------------
# Resolve directories
# --------------------------------------------------
if args.base_dir:
    BASE = args.base_dir
    PRES_DIR = args.presence_dir or os.path.join(BASE, "genepresenceabsence")
    METRIC_DIR = args.metric_dir or os.path.join(BASE, "metrics")
else:
    PRES_DIR = args.presence_dir
    METRIC_DIR = args.metric_dir

OUT_TEX = args.out_tex

print("[INFO] Generating LaTeX report")

# --------------------------------------------------
# Find genes
# --------------------------------------------------
presence_files = glob.glob(os.path.join(PRES_DIR, "*_gene_presence.png"))

genes = sorted([
    os.path.basename(f).replace("_gene_presence.png", "")
    for f in presence_files
])

print("[INFO] Genes found:", len(genes))

# --------------------------------------------------
# Escape LaTeX characters
# --------------------------------------------------
def esc(text):
    return text.replace("_", r"\_")

# --------------------------------------------------
# Write LaTeX
# --------------------------------------------------
with open(OUT_TEX, "w") as tex:

    tex.write(r"""
\documentclass[a4paper]{article}
\usepackage{graphicx}
\usepackage{longtable}
\usepackage[margin=1in]{geometry}
\usepackage{hyperref}

\begin{document}

\begin{center}
{\LARGE \textbf{Synteny Verification Report}}\\
\vspace{0.5cm}
{\large Reference species: \textit{Homo sapiens}}
\end{center}

\newpage
""")

    # --------------------------------------------------
    # Index
    # --------------------------------------------------
    tex.write(r"\section*{Index}" + "\n")

    for g in genes:
        tex.write(f"\\hyperref[sec:{g}]{{{esc(g)}}}\\\\\n")

    tex.write("\n\\newpage\n")

    # --------------------------------------------------
    # Per-gene section
    # --------------------------------------------------
    for gene in genes:

        pres_plot = os.path.join(
            PRES_DIR,
            f"{gene}_gene_presence.png"
        )

        metric_file = os.path.join(
            METRIC_DIR,
            f"{gene}_metrics.tsv"
        )

        tex.write(f"\\section*{{{esc(gene)}}}\n")
        tex.write(f"\\label{{sec:{gene}}}\n")

        # Presence plot
        if os.path.exists(pres_plot):

            tex.write(r"""
\subsection*{Gene Presence}
\begin{center}
""")

            tex.write(
                f"\\includegraphics[width=\\textwidth]{{{pres_plot}}}\n"
            )

            tex.write(r"\end{center}")

        tex.write("\n")

        # Metrics table
        if os.path.exists(metric_file):

            df = pd.read_csv(metric_file, sep="\t")

            tex.write(r"""
\subsection*{Synteny Metrics}
\begin{center}
\begin{tabular}{lccc}
Species & Presence & Order & Score \\
\hline
""")

            for _, row in df.iterrows():

                sp = esc(str(row["Species"]))
                p = row["Presence_Score"]
                o = row["Order_Score"]
                s = row["Synteny_Score"]

                tex.write(
                    f"{sp} & {p} & {o} & {s} \\\\\n"
                )

            tex.write(r"""
\end{tabular}
\end{center}
""")

        tex.write("\n\\newpage\n")

    tex.write(r"\end{document}")

print("[OK] LaTeX file created:", OUT_TEX)