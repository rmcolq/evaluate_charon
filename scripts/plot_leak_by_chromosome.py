#!/usr/bin/env python3
"""How leaked reads distribute over chromosomes, for each classifier.

Input: leak_reads.tsv from leak_regions.py; samples are pooled (run once per read type). Top panel: leaked reads
per chromosome (grouped bars, one per classifier). Bottom panel: the same per 100 Mb of chromosome, so the
long chromosomes do not dominate; a uniform spread over the genome would give a flat line. Writes the
plot and `leak_by_chromosome.tsv` (counts, per-100Mb rates and a chi-square test against length-proportional
expectation for each classifier).

  scripts/plot_leak_by_chromosome.py --reads results/leak_regions/leak_reads.tsv --samples HG02583 \
      --classifiers charon deacon deacon_sensitive nohuman --ref_bed resources/T2T_ref.bed -o by_chromosome.png
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chisquare

MITO = "JAGYVL020000058.1"  # human mitochondrion contig in the default refs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--classifiers", nargs="+", required=True)
    ap.add_argument("--ref_bed", default="resources/T2T_ref.bed")
    ap.add_argument("--title", default=None)
    ap.add_argument("-o", "--output", default="leak_by_chromosome.png")
    a = ap.parse_args()

    names, lengths = [], {}
    with open(a.ref_bed) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and (p[0].startswith("NC_0609") or p[0] == MITO):
                lengths[p[0]] = int(p[2])
                desc = p[3]
                lab = ("chrM" if p[0] == MITO else "chr" + desc.split("chromosome ")[1].split(",")[0])
                names.append((p[0], lab))
    refs = [r for r, _ in names]
    labels = [l for _, l in names]

    reads = pd.read_csv(a.reads, sep="\t")
    reads = reads[reads["sample"].isin(a.samples) & (reads["destination"] == "T2T")]
    reads = reads.assign(ref=reads["ref"].replace({"NC_012920.1": MITO, "CP068254.1": MITO}))

    rows = []
    counts = {}
    for clf in a.classifiers:
        c = reads[reads["classifier"] == clf]["ref"].value_counts()
        counts[clf] = np.array([int(c.get(r, 0)) for r in refs])
    nuclear = [i for i, r in enumerate(refs) if r != MITO]  # test on nuclear chromosomes only
    length_arr = np.array([lengths[r] for r in refs], dtype=float)
    for clf in a.classifiers:
        obs = counts[clf]
        n = int(obs.sum())
        p_chi = float("nan")
        n_nuc = int(obs[nuclear].sum())
        if n_nuc >= 20:
            exp = n_nuc * length_arr[nuclear] / length_arr[nuclear].sum()
            p_chi = chisquare(obs[nuclear], exp).pvalue
        for r, l, k, ln in zip(refs, labels, obs, length_arr):
            rows.append({"classifier": clf, "ref": r, "chromosome": l, "n_reads": int(k),
                         "reads_per_100Mb": k / ln * 1e8, "chisq_p_vs_length": p_chi})
    pd.DataFrame(rows).to_csv(os.path.splitext(a.output)[0] + ".tsv", sep="\t", index=False)

    colours = dict(zip(a.classifiers, plt.cm.Dark2(np.arange(len(a.classifiers)) % 8)))
    x = np.arange(len(refs))
    w = 0.8 / len(a.classifiers)
    fig, axes = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
    for j, clf in enumerate(a.classifiers):
        off = (j - (len(a.classifiers) - 1) / 2) * w
        total = int(counts[clf].sum())
        axes[0].bar(x + off, counts[clf], w, color=colours[clf], label=f"{clf} (n={total})")
        axes[1].bar(x + off, counts[clf] / length_arr * 1e8, w, color=colours[clf])
    axes[0].set_ylabel("leaked reads")
    axes[1].set_ylabel("leaked reads per 100 Mb")
    axes[1].set_yscale("symlog", linthresh=1)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=60, fontsize=8)
    axes[0].legend(fontsize=8, frameon=False)
    axes[1].text(0.995, 0.97, "chrM is a 16.5 kb contig, so its per-100Mb rate is large by construction", transform=axes[1].transAxes,
                 ha="right", va="top", fontsize=7, color="#666666")
    axes[0].set_title(f"{a.title + ': ' if a.title else ''}leaked reads per chromosome (pooled: {', '.join(a.samples)})",
                      fontsize=10, loc="left")
    for ax in axes:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.output, dpi=150)
    print("wrote", a.output)


if __name__ == "__main__":
    main()
