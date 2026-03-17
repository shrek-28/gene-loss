#!/usr/bin/env python3
import os
import argparse
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description="Generate LaTeX PDF for Klumpy Alignment Verification")
    parser.add_argument("--plots_dir", required=True, help="Directory containing plots subfolders per gene")
    parser.add_argument("--metrics_dir", required=True, help="Directory containing metrics subfolders per gene")
    parser.add_argument("--flanking_tsv", required=True, help="Path to flanking_regions.tsv")
    parser.add_argument("--output_tex", required=True, help="Output LaTeX file path")
    parser.add_argument("--organism_name", required=True, help="Organism name to appear in title page")
    parser.add_argument("--output_pdf", required=True, help="Output PDF file path")
    return parser.parse_args()

def escape_latex(s):
    """Escape LaTeX special characters"""
    import re
    s = s.replace("_", " ")
    s = s.replace("&", "\\&")
    s = s.replace("%", "\\%")
    s = s.replace("$", "\\$")
    s = s.replace("#", "\\#")
    s = s.replace("{", "\\{")
    s = s.replace("}", "\\}")
    s = s.replace("~", "\\textasciitilde{}")
    s = s.replace("^", "\\textasciicircum{}")
    s = s.replace("\\", "\\textbackslash{}")
    s = re.sub(r"\s+", " ", s)
    return s.strip()

def tsv_to_latex_table(tsv_file, table_type="common", transpose=False):
    """Convert TSV to LaTeX table that never exceeds page width"""

    if not os.path.exists(tsv_file):
        return f"% {tsv_file} not found\n"

    with open(tsv_file) as f:
        lines = [line.strip() for line in f if line.strip()]
    if not lines:
        return f"% {tsv_file} empty\n"

    headers = [escape_latex(x) for x in lines[0].split("\t")]
    rows = [[escape_latex(x) for x in line.split("\t")] for line in lines[1:]]

    if transpose:
        transposed = []
        for i in range(len(headers)):
            transposed.append([headers[i]] + [row[i] for row in rows])
        headers = ["Metric"] + [f"Sample {i+1}" for i in range(len(rows))]
        rows = transposed

    col_format = " | ".join(["l"] * len(headers))

    latex = "\\begin{table}[H]\n"
    latex += "\\centering\n"
    latex += f"\\caption{{{table_type.capitalize()} metrics}}\n"

    latex += "\\begin{adjustbox}{max width=\\textwidth}\n"
    latex += f"\\begin{{tabular}}{{{col_format}}}\n"
    latex += "\\hline\n"

    latex += " & ".join(headers) + " \\\\\n"
    latex += "\\hline\n"

    for row in rows:
        latex += " & ".join(row) + " \\\\\n"

    latex += "\\hline\n"
    latex += "\\end{tabular}\n"
    latex += "\\end{adjustbox}\n"
    latex += "\\end{table}\n\n"

    return latex

def plot_sort_key(plot_file):
    """Sort plots descending by flank size: 500 > 250 > 100 > 50"""
    for size in ["500", "250", "100", "50"]:
        if f"_flank_{size}" in plot_file:
            return -int(size)  # negative to sort descending
    return 0  # default goes last

def generate_gene_section(gene, args):
    """Return LaTeX code string for one focal gene"""
    gene_plot_dir = os.path.join(args.plots_dir, gene)
    gene_metric_dir = os.path.join(args.metrics_dir, gene)
    latex = f"\\section{{{escape_latex(gene)}}}\n"

    # Plots sorted by flank size descending
    plot_files = [f for f in os.listdir(gene_plot_dir) if "_alignment_plot.pdf" in f]
    plot_files.sort(key=plot_sort_key)
    for plot_file in plot_files:
        if "_flank_" in plot_file:
            threshold = plot_file.split("_flank_")[1].split("_alignment_plot.pdf")[0] + " kb"
        else:
            threshold = "default"
        plot_path = os.path.join(gene_plot_dir, plot_file).replace("\\","/")
        latex += f"""
\\begin{{figure}}[H]
\\centering
\\includegraphics[width=0.9\\textwidth,height=0.8\\textheight,keepaspectratio]{{{plot_path}}}
\\caption{{{escape_latex(gene)}: {threshold} threshold}}
\\end{{figure}}
\\vspace{{1em}}
"""

    # Tables: common and threshold
    common_metrics_file = os.path.join(gene_metric_dir, f"{gene}_common_metrics.tsv")
    threshold_metrics_file = os.path.join(gene_metric_dir, f"{gene}_threshold_metrics.tsv")

    latex += "\\subsection*{Common Metrics}\n"
    latex += tsv_to_latex_table(common_metrics_file, table_type="common", transpose=False)
    latex += "\n\\subsection*{Threshold Metrics}\n"
    latex += tsv_to_latex_table(threshold_metrics_file, table_type="threshold", transpose=True)
    latex += "\n\\newpage\n"

    return latex

def generate_appendix(flanking_tsv):
    latex = "\\appendix\n\\section{Flanking Regions Table}\n"
    latex += tsv_to_latex_table(flanking_tsv, table_type="flanking regions")
    return latex

def main():
    args = parse_args()

    genes = [d for d in os.listdir(args.plots_dir) if os.path.isdir(os.path.join(args.plots_dir, d))]
    genes.sort()

    latex_content = f"""
\\documentclass[a4paper,12pt]{{article}}
\\usepackage[margin=1cm]{{geometry}}
\\usepackage{{graphicx}}
\\usepackage{{caption}}
\\usepackage{{multicol}}
\\usepackage{{hyperref}}
\\usepackage{{float}}
\\usepackage{{fancyhdr}}
\\usepackage{{longtable}}
\\usepackage{{adjustbox}}

\\pagestyle{{fancy}}
\\fancyhf{{}}
\\rhead{{\\thepage}}

\\title{{Klumpy Alignment Verification Results}}
\\author{{{escape_latex(args.organism_name)}}}

\\begin{{document}}
\\maketitle
\\newpage
\\tableofcontents
\\newpage
"""

    for gene in genes:
        latex_content += generate_gene_section(gene, args)

    # latex_content += generate_appendix(args.flanking_tsv)
    latex_content += "\n\\end{document}"

    with open(args.output_tex, "w") as f:
        f.write(latex_content)

    # Compile PDF
    subprocess.run(["pdflatex", "-interaction=nonstopmode", args.output_tex])
    subprocess.run(["pdflatex", "-interaction=nonstopmode", args.output_tex])
    print(f"PDF generated: {args.output_pdf}")

if __name__ == "__main__":
    main()