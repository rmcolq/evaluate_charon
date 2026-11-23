#!/usr/bin/env python3
"""Where do microbial-called reads that map to the human (T2T) reference come from?

Reads the `<sample>.<classifier>.microbial.sam` files published in `<outdir>/<sample>/mapped/` (only reads that
the classifier called microbial are in them) and writes, per sample and classifier:

  * leak_reads.tsv   one row per leaked read: locus, MAPQ, divergence, aligned length, query length
  * leak_loci.tsv    reads merged into loci (alignments within --merge_dist bp), with how many reads support each
  * leak_bins.tsv    reads per fixed-size genome bin (--bin_size), for an ideogram-style plot
  * leak_summary.tsv counts by destination (T2T, EBV, other contaminant, unmapped) per sample and classifier

Standard library only. If matplotlib is available, an ideogram-style figure is written too.
Optionally restrict to reads another classifier also processed (--common_ids), because Charon silently
drops reads below its min quality / length, which the other tools keep.

  scripts/leak_regions.py --outdir out --samples HG02583 HG02922 --classifiers charon deacon nohuman
"""
import argparse
import csv
import os
import re
from collections import defaultdict

EBV = "NC_007605.1"
T2T_PREFIX = "NC_0609"  # CHM13v2.0 RefSeq accessions NC_060925-NC_060948
# chrM. The default refs carry the human mitochondrion as the HPRC contig JAGYVL020000058.1 (16,544 bp);
# the CHM13 and rCRS names are accepted in case the refs are changed.
MITO = ("JAGYVL020000058.1", "CP068254.1", "NC_012920.1")

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_sam(path):
    """Yield (read_id, ref, pos, end, mapq, divergence, aligned_len, query_len, is_primary) for each line."""
    ref_lengths = {}
    with open(path) as f:
        for line in f:
            if line.startswith("@SQ"):
                fields = dict(x.split(":", 1) for x in line.rstrip().split("\t")[1:])
                ref_lengths[fields["SN"]] = int(fields["LN"])
                continue
            if line.startswith("@"):
                continue
            p = line.rstrip("\n").split("\t")
            flag, ref, pos, mapq, cigar = int(p[1]), p[2], int(p[3]), int(p[4]), p[5]
            if ref == "*":
                continue
            ref_len = qry_len = aln_len = 0
            for n, op in CIGAR_RE.findall(cigar):
                n = int(n)
                if op in "M=X":
                    ref_len += n; aln_len += n; qry_len += n
                elif op in "DN":
                    ref_len += n
                elif op in "IS":
                    qry_len += n
                elif op == "H":
                    qry_len += n
            de = next((float(x[5:]) for x in p[11:] if x.startswith("de:f:")), None)
            yield p[0], ref, pos, pos + ref_len - 1, mapq, de, aln_len, qry_len, not (flag & 0x900), ref_lengths
    return


def destination(ref):
    if ref == EBV:
        return "EBV"
    if ref.startswith(T2T_PREFIX) or ref in MITO:
        return "T2T"
    return "other_contaminant"


def cluster(rows, dist):
    """Merge alignments (ref, start, end, read_id) into loci; rows are sorted by (ref, start)."""
    loci = []
    for ref, start, end, rid in sorted(rows):
        if loci and loci[-1]["ref"] == ref and start <= loci[-1]["end"] + dist:
            loci[-1]["end"] = max(loci[-1]["end"], end)
            loci[-1]["reads"].add(rid)
        else:
            loci.append({"ref": ref, "start": start, "end": end, "reads": {rid}})
    return loci


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="out", help="pipeline output directory (contains <sample>/mapped/)")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--classifiers", nargs="+", required=True)
    ap.add_argument("--results", default="results/leak_regions", help="where to write tables")
    ap.add_argument("--merge_dist", type=int, default=1000, help="merge alignments closer than this (bp)")
    ap.add_argument("--bin_size", type=int, default=1_000_000)
    ap.add_argument("--common_ids", default=None,
                    help="TSV/CSV/.out file whose second column holds read ids; keep only those reads (e.g. the reads "
                         "Charon processed, `<sample>_charon.out`). Applied per sample: use {sample} as a placeholder")
    args = ap.parse_args()
    os.makedirs(args.results, exist_ok=True)

    reads_out = open(os.path.join(args.results, "leak_reads.tsv"), "w", newline="")
    loci_out = open(os.path.join(args.results, "leak_loci.tsv"), "w", newline="")
    bins_out = open(os.path.join(args.results, "leak_bins.tsv"), "w", newline="")
    sum_out = open(os.path.join(args.results, "leak_summary.tsv"), "w", newline="")
    rw, lw, bw, sw = (csv.writer(x, delimiter="\t") for x in (reads_out, loci_out, bins_out, sum_out))
    rw.writerow(["sample", "classifier", "read_id", "destination", "ref", "start", "end", "mapq",
                 "divergence", "aligned_len", "query_len", "aligned_frac"])
    lw.writerow(["sample", "classifier", "ref", "start", "end", "span", "n_reads"])
    bw.writerow(["sample", "classifier", "ref", "ref_length", "bin_start", "n_reads"])
    sw.writerow(["sample", "classifier", "microbial_reads_mapped", "T2T", "EBV", "other_contaminant"])

    for sample in args.samples:
        keep = None
        if args.common_ids:
            keep = set()
            with open(args.common_ids.format(sample=sample)) as f:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) > 1:
                        keep.add(p[1])
        for clf in args.classifiers:
            path = os.path.join(args.outdir, sample, "mapped", f"{sample}.{clf}.microbial.sam")
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                print(f"skip {path}: missing or empty (was --evaluate_microbial on?)")
                continue
            best = {}  # read_id -> best (highest mapq, then longest aligned) primary alignment
            ref_lengths = {}
            for rid, ref, pos, end, mapq, de, aln, qlen, primary, ref_lengths in parse_sam(path):
                if keep is not None and rid not in keep:
                    continue
                if not primary:
                    continue
                cand = (mapq, aln, ref, pos, end, de, qlen)
                if rid not in best or cand[:2] > best[rid][:2]:
                    best[rid] = cand
            counts = defaultdict(int)
            t2t_rows = []
            bins = defaultdict(int)
            for rid, (mapq, aln, ref, pos, end, de, qlen) in best.items():
                dest = destination(ref)
                counts[dest] += 1
                rw.writerow([sample, clf, rid, dest, ref, pos, end, mapq, "" if de is None else de, aln, qlen,
                             round(aln / qlen, 3) if qlen else ""])
                if dest == "T2T":
                    t2t_rows.append((ref, pos, end, rid))
                    bins[(ref, ((pos + end) // 2) // args.bin_size * args.bin_size)] += 1
            for l in cluster(t2t_rows, args.merge_dist):
                lw.writerow([sample, clf, l["ref"], l["start"], l["end"], l["end"] - l["start"] + 1, len(l["reads"])])
            for (ref, b), n in sorted(bins.items()):
                bw.writerow([sample, clf, ref, ref_lengths.get(ref, ""), b, n])
            sw.writerow([sample, clf, len(best), counts["T2T"], counts["EBV"], counts["other_contaminant"]])
            print(f"{sample} {clf}: {len(best)} mapped reads, T2T={counts['T2T']} EBV={counts['EBV']} "
                  f"other={counts['other_contaminant']}")
    for x in (reads_out, loci_out, bins_out, sum_out):
        x.close()


if __name__ == "__main__":
    main()
