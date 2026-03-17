#!/usr/bin/env python3

import argparse
import subprocess
import os
import sys

def main():

    parser = argparse.ArgumentParser(
        description="Run klumpy scan_alignments for a single species."
    )

    # Required inputs
    parser.add_argument("--species", required=True,
                        help="Species name (used for output naming)")

    parser.add_argument("--bam", required=True,
                        help="Sorted and indexed BAM file")

    parser.add_argument("--outdir", required=True,
                        help="Output directory")

    # Optional inputs
    parser.add_argument("--annotation", default=None,
                        help="GTF/GFF annotation file")

    parser.add_argument("--threads", type=int, default=1)

    parser.add_argument("--min_len", type=int, default=2000)
    parser.add_argument("--min_percent", type=float, default=50)

    parser.add_argument("--num_of_groups", type=int, default=3)

    parser.add_argument("--window_size", type=int, default=50000)
    parser.add_argument("--window_step", type=int, default=25000)

    parser.add_argument("--limit", type=int, default=10000)

    parser.add_argument("--flag_excess_groups", action="store_true")

    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    output_prefix = os.path.join(args.outdir, args.species)

    cmd = [
        "klumpy", "scan_alignments",
        "--alignment_map", args.bam,
        "--threads", str(args.threads),
        "--min_len", str(args.min_len),
        "--min_percent", str(args.min_percent),
        "--num_of_groups", str(args.num_of_groups),
        "--window_size", str(args.window_size),
        "--window_step", str(args.window_step),
        "--limit", str(args.limit),
    ]

    if args.annotation:
        cmd.extend(["--annotation", args.annotation])

    if args.flag_excess_groups:
        cmd.append("--flag_excess_groups")

    print("Running command:")
    print(" ".join(cmd))
    print()

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError:
        print("Error: scan_alignments failed", file=sys.stderr)
        sys.exit(1)

    print("Finished scan_alignments for:", args.species)
    print("Outputs written to:", args.outdir)


if __name__ == "__main__":
    main()