#!/usr/bin/env python
"""Write one row per read of a SAM file: where it maps, how much of the read is aligned, and whether it maps in parts.
Uses the same columns as the full tables of evaluate.py, so reads which were never classified (e.g. charon's
unclassified reads) can be compared with the ones which were."""
import argparse
import sys

from evaluate import load_map_info_from_sam


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sam", required=True, help="SAM file from minimap2")
    parser.add_argument("--out", required=True, help="output CSV (gzipped if it ends in .gz)")
    args = parser.parse_args()
    df = load_map_info_from_sam(args.sam)
    df.to_csv(args.out, index=False)
    sys.stderr.write(f"Wrote {df.shape[0]} reads to {args.out}\n")


if __name__ == "__main__":
    main()
