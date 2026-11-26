#!/usr/bin/env python3
"""Per-chromosome plots of where leaked reads map on the T2T reference, and what sequence is there.

Input: leak_reads.tsv from leak_regions.py. Reads from all `--samples` are pooled, so run it once per read type
(nanopore, illumina). Writes one PNG per chromosome (chr1-22, X, Y, M) into `--outdir`, plus all of them in one
PDF. Each figure has one copy of the annotation tracks on top, drawn on the same bins as the read heatmaps:

  * repeat content: stacked bars, the share of each bin that is LINE / SINE / LTR / DNA / satellite /
    simple + low complexity / rRNA / small RNA / other repeat (white = no repeat, i.e. unique sequence)
  * genes: share of each bin that is genic (transcript span) and exonic (RefSeq)
  * CenSat + rDNA, and segmental duplications, as thin strips at full resolution

Below them, one heatmap track per classifier: bins coloured by reads per bin (colour gets stronger with more reads,
one scale shared by all classifiers and chromosomes). A classifier with no leaked reads still gets its (empty) track.

Also writes `enrichment.tsv`: for each classifier and annotation, the share of leaked-read aligned bases that fall
in the annotation versus the share of the genome it covers (enrichment > 1: leaks are over-represented there).

  scripts/plot_leak_ideogram.py --reads results/leak_regions/leak_reads.tsv \
      --samples HG02583 HG02922 NA18983 --classifiers charon deacon deacon_sensitive nohuman \
      --annotations resources/annotations --outdir results/ideogram_nanopore

The annotation files come from the UCSC hub for T2T-CHM13v2.0. Any that are missing from `--annotations` are
downloaded there on first use (about 215 MB) and checked against UCSC's md5 checksums; `--no_download` turns this off.
"""
import argparse
import gzip
import hashlib
import os
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

MITO = "JAGYVL020000058.1"  # human mitochondrion contig in the default refs (labelled chrM)

# Strip annotations drawn at full resolution (CenSat classes, rDNA, segdups)
STRIP_COLOURS = {
    "alpha satellite": "#d62728", "HSat1": "#2ca02c", "HSat2": "#1f77b4", "HSat3": "#17becf",
    "beta satellite": "#e377c2", "gamma satellite": "#9467bd", "rDNA": "#000000", "segdup": "#ff7f0e",
}
# Repeat classes (RepeatMasker class collapsed) drawn as stacked bars per bin
REPEAT_COLOURS = {
    "LINE": "#4c78a8", "SINE": "#f58518", "LTR": "#54a24b", "DNA": "#b279a2", "satellite": "#e45756",
    "simple/low complexity": "#9d755d", "rRNA": "#000000", "small RNA": "#bab0ac", "other repeat": "#79706e",
}


UCSC_HUB = "https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4"
# (file, directory under UCSC_HUB, md5 from UCSC_HUB/md5sum.txt)
ANNOTATION_FILES = [
    ("GCA_009914755.4.chromAlias.txt", "", "3c7f3f0529a46d9da6fee6aee70a1fed"),
    ("t2t_censat_CHM13v2.0_trackv2.0.10col.bed", "bbi/GCA_009914755.4_T2T-CHM13v2.0.censat",
     "b50c88e0e417ba352f7b8a230e5057be"),
    ("chm13v1.1.rdna_model.bed", "bbi/GCA_009914755.4_T2T-CHM13v2.0.rdnaModel", "057746b63d5a116aee4e016deaf53dd5"),
    ("T2T-CHM13v2.SDs.bed.gz", "bbi/GCA_009914755.4_T2T-CHM13v2.0.sedefSegDups", "5d3f6592942fd76a8acbaa177e0993cb"),
    ("chm13v2.0_rmsk.out.gz", "bbi/GCA_009914755.4_T2T-CHM13v2.0.t2tRepeatMasker", "eac086f914d584e9cdc190e42747d377"),
    ("GCF_009914755.1_T2T-CHM13v2.0.110.20220412.ncbiRefSeq.gtf.gz", "genes", "d5fb63080182480fa9b8e13f27fcaa34"),
]


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_annotations(ann_dir):
    """Download any missing annotation file into ann_dir, checking its md5 before keeping it."""
    os.makedirs(ann_dir, exist_ok=True)
    for name, subdir, md5 in ANNOTATION_FILES:
        dest = os.path.join(ann_dir, name)
        if os.path.exists(dest):
            continue
        url = "/".join(x for x in (UCSC_HUB, subdir, name) if x)
        print(f"downloading {url} ...", flush=True)
        tmp = dest + ".part"
        urllib.request.urlretrieve(url, tmp)
        got = md5sum(tmp)
        if got != md5:
            os.remove(tmp)
            raise SystemExit(f"md5 mismatch for {url}: expected {md5}, got {got}")
        os.replace(tmp, dest)


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def censat_class(name):
    n = name.lower()
    if n.startswith(("hor", "dhor", "mon")):
        return "alpha satellite"
    for k, v in (("hsat1", "HSat1"), ("hsat2", "HSat2"), ("hsat3", "HSat3"), ("bsat", "beta satellite"),
                 ("gsat", "gamma satellite"), ("rdna", "rDNA")):
        if n.startswith(k):
            return v
    return None  # ct_ (transition) and generic censat_ regions are not drawn


def repeat_class(cls):
    c = cls.split("/")[0]
    if c in ("Simple_repeat", "Low_complexity"):
        return "simple/low complexity"
    if c in ("LINE", "SINE", "LTR", "DNA", "rRNA"):
        return c
    if c == "Satellite":
        return "satellite"
    if c in ("tRNA", "snRNA", "scRNA", "srpRNA"):
        return "small RNA"
    return "other repeat"


def load_alias(path):
    """Return {ucsc name: refseq-or-genbank accession used in the SAMs}, {genbank: ucsc name}."""
    to_acc, gb_to_ucsc = {}, {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 5:
                to_acc[p[4]] = p[1] or p[0]
                gb_to_ucsc[p[0]] = p[4]
    return to_acc, gb_to_ucsc


def merge(starts, ends):
    """Merge overlapping 0-based half-open intervals; returns sorted (starts, ends) arrays."""
    if len(starts) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    o = np.argsort(starts, kind="stable")
    s, e = starts[o].astype(np.int64), ends[o].astype(np.int64)
    run_end = np.maximum.accumulate(e)
    new = np.ones(len(s), dtype=bool)
    new[1:] = s[1:] > run_end[:-1]
    idx = np.flatnonzero(new)
    ms = s[idx]
    me = np.maximum.reduceat(e, idx)
    return ms, me


def to_layers(df, chrom_col, start_col, end_col, layer_col, to_acc):
    """DataFrame of intervals -> {layer: {accession: (starts, ends)}} (merged)."""
    out = {}
    for (layer, chrom), g in df.groupby([layer_col, chrom_col]):
        acc = to_acc.get(chrom)
        if acc:
            out.setdefault(layer, {})[acc] = merge(g[start_col].to_numpy(), g[end_col].to_numpy())
    return out


def covered(ivs, x):
    """Bases of the merged intervals ivs=(starts, ends) lying in [0, x), vectorised over x."""
    s, e = ivs
    x = np.asarray(x, dtype=np.int64)
    ps = np.concatenate([[0], np.cumsum(s)])
    pe = np.concatenate([[0], np.cumsum(e)])
    i = np.searchsorted(s, x, side="left")
    j = np.searchsorted(e, x, side="left")
    return (i * x - ps[i]) - (j * x - pe[j])


def overlap(ivs, lo, hi):
    return covered(ivs, hi) - covered(ivs, lo)


def load_annotations(ann_dir, to_acc, gb_to_ucsc):
    """Return (strips, repeats, genes): each {layer: {accession: (starts, ends)}}."""
    strips, repeats, genes = {}, {}, {}
    p = os.path.join(ann_dir, "t2t_censat_CHM13v2.0_trackv2.0.10col.bed")
    if os.path.exists(p):
        df = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["c", "s", "e", "n"],
                         skiprows=1)
        df["layer"] = df["n"].map(censat_class)
        strips.update(to_layers(df.dropna(subset=["layer"]), "c", "s", "e", "layer", to_acc))
    p = os.path.join(ann_dir, "chm13v1.1.rdna_model.bed")
    if os.path.exists(p):
        df = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2], names=["c", "s", "e"])
        df["layer"] = "rDNA"
        for k, v in to_layers(df, "c", "s", "e", "layer", to_acc).items():
            strips.setdefault(k, {}).update(v)
    p = os.path.join(ann_dir, "T2T-CHM13v2.SDs.bed.gz")
    if os.path.exists(p):
        df = pd.read_csv(p, sep="\t", header=None, skiprows=1, usecols=[0, 1, 2], names=["c", "s", "e"])
        df["layer"] = "segdup"
        strips.update(to_layers(df, "c", "s", "e", "layer", to_acc))
    p = os.path.join(ann_dir, "chm13v2.0_rmsk.out.gz")
    if os.path.exists(p):
        print("reading RepeatMasker ...", flush=True)
        df = pd.read_csv(p, sep=r"\s+", skiprows=3, header=None, usecols=[4, 5, 6, 10],
                         names=["c", "s", "e", "cls"], dtype={"c": str, "cls": str})
        df["s"] = df["s"] - 1
        df["layer"] = df["cls"].map(repeat_class)
        repeats.update(to_layers(df, "c", "s", "e", "layer", to_acc))
        df["layer"] = "any repeat"
        repeats.update(to_layers(df, "c", "s", "e", "layer", to_acc))
    p = os.path.join(ann_dir, "GCF_009914755.1_T2T-CHM13v2.0.110.20220412.ncbiRefSeq.gtf.gz")
    if os.path.exists(p):
        print("reading genes ...", flush=True)
        df = pd.read_csv(p, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4], names=["c", "f", "s", "e"])
        df["c"] = df["c"].map(gb_to_ucsc)
        df = df.dropna(subset=["c"])
        df["s"] = df["s"] - 1
        for feat, layer in (("transcript", "genic"), ("exon", "exonic")):
            sub = df[df["f"] == feat].copy()
            sub["layer"] = layer
            genes.update(to_layers(sub, "c", "s", "e", "layer", to_acc))
    return strips, repeats, genes


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reads", required=True, help="leak_reads.tsv from leak_regions.py")
    ap.add_argument("--samples", nargs="+", required=True, help="samples to pool (one read type per run)")
    ap.add_argument("--classifiers", nargs="+", required=True, help="tracks, in order (empty if no leaks)")
    ap.add_argument("--annotations", default="resources/annotations",
                    help="directory of the UCSC annotation files; missing ones are downloaded into it")
    ap.add_argument("--no_download", action="store_true",
                    help="use only the annotation files already in --annotations")
    ap.add_argument("--ref_bed", default="resources/T2T_ref.bed", help="chromosome lengths (name, 0, length)")
    ap.add_argument("--bin_size", type=int, default=1_000_000,
                    help="bin for the heatmaps and bar tracks (chromosomes under 5 Mb use 50 bins)")
    ap.add_argument("--title", default=None, help="prefix for each figure title, e.g. 'nanopore'")
    ap.add_argument("--outdir", default="results/ideogram")
    ap.add_argument("--chromosomes", nargs="*", default=None, help="e.g. chr1 chrX chrM (default: all)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if not a.no_download:
        fetch_annotations(a.annotations)

    to_acc, gb_to_ucsc = load_alias(os.path.join(a.annotations, "GCA_009914755.4.chromAlias.txt"))
    rev = {v: k for k, v in to_acc.items()}
    rev[MITO] = "chrM"
    strips, repeats, genes = load_annotations(a.annotations, to_acc, gb_to_ucsc)
    have_rep = any(k in repeats for k in REPEAT_COLOURS)
    have_genes = "genic" in genes

    reads = pd.read_csv(a.reads, sep="\t")
    reads = reads[reads["sample"].isin(a.samples) & (reads["destination"] == "T2T")].copy()
    reads["ref"] = reads["ref"].replace({"NC_012920.1": MITO, "CP068254.1": MITO})
    reads["mid"] = (reads["start"] + reads["end"]) / 2

    lengths = {}
    if os.path.exists(a.ref_bed):
        with open(a.ref_bed) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) >= 3:
                    lengths[p[0]] = int(p[2])
    else:
        print(f"warning: {a.ref_bed} not found; chromosome lengths fall back to the last annotated base")
    lengths.setdefault(MITO, 16544)
    chroms = [c for c in rev if c.startswith("NC_0609") or c == MITO]
    chroms = sorted(chroms, key=lambda c: (0, int(rev[c][3:])) if rev[c][3:].isdigit() else (1, "XYM".index(rev[c][3:])))
    if a.chromosomes:
        chroms = [c for c in chroms if rev[c] in a.chromosomes]

    def length_of(c):
        if c in lengths:
            return lengths[c]
        return max((int(v[c][1].max()) for lay in (strips, repeats, genes) for v in lay.values() if c in v), default=1)

    def bin_width(ln):
        return a.bin_size if ln > 5_000_000 else max(1, ln // 50)

    # one colour scale for every chromosome and classifier, so the same colour always means the same read count
    gmax = 1
    for c in chroms:
        ln = length_of(c)
        edges = np.arange(0, ln + bin_width(ln), bin_width(ln))
        for clf in a.classifiers:
            mids = reads[(reads["classifier"] == clf) & (reads["ref"] == c)]["mid"]
            if len(mids):
                gmax = max(gmax, int(np.histogram(mids, bins=edges)[0].max()))
    norm = mcolors.Normalize(vmin=0, vmax=gmax)
    cmap = matplotlib.colormaps["viridis"]

    pdf = PdfPages(os.path.join(a.outdir, "leak_chromosomes.pdf"))
    for c in chroms:
        ln = length_of(c)
        scale, unit = (1e6, "Mb") if ln > 1_000_000 else (1e3, "kb")
        bin_bp = bin_width(ln)
        edges = np.arange(0, ln + bin_bp, bin_bp)
        centers = (edges[:-1] + bin_bp / 2) / scale
        widths = np.diff(np.minimum(edges, ln))  # last bin may be partial
        n_t = len(a.classifiers)
        n_top = 2 + (1 if have_rep else 0) + (1 if have_genes else 0)
        fig, axes = plt.subplots(n_top + n_t, 1, figsize=(12, 1.2 + 0.8 * (n_top + n_t)), sharex=True,
                                 gridspec_kw={"height_ratios": ([1.6] if have_rep else []) + ([1.2] if have_genes else [])
                                                               + [0.6, 0.6] + [1] * n_t})
        ax_i = 0
        legend_items = []
        if have_rep:
            ax = axes[ax_i]; ax_i += 1
            bottom = np.zeros(len(widths))
            for cls, col in REPEAT_COLOURS.items():
                if cls in repeats and c in repeats[cls]:
                    frac = overlap(repeats[cls][c], edges[:-1], np.minimum(edges[1:], ln)) / widths
                    ax.bar(centers, frac, width=bin_bp / scale, bottom=bottom, color=col, lw=0)
                    bottom += frac
                    legend_items.append((cls, col))
            ax.set_ylim(0, 1)
            ax.set_ylabel("repeat content\n(share of bin)", rotation=0, ha="right", va="center", fontsize=7)
        if have_genes:
            ax = axes[ax_i]; ax_i += 1
            if c in genes.get("genic", {}):
                g = overlap(genes["genic"][c], edges[:-1], np.minimum(edges[1:], ln)) / widths
                ax.bar(centers, g, width=bin_bp / scale, color="#9ecae1", lw=0)
                legend_items.append(("genic", "#9ecae1"))
            if c in genes.get("exonic", {}):
                x = overlap(genes["exonic"][c], edges[:-1], np.minimum(edges[1:], ln)) / widths
                ax.bar(centers, x, width=bin_bp / scale, color="#08519c", lw=0)
                legend_items.append(("exonic", "#08519c"))
            ax.set_ylim(0, 1)
            ax.set_ylabel("genes\n(share of bin)", rotation=0, ha="right", va="center", fontsize=7)
        sat, sd = axes[ax_i], axes[ax_i + 1]
        ax_i += 2
        for cls, per_chr in strips.items():
            ax = sd if cls == "segdup" else sat
            if c in per_chr:
                for s_, e_ in zip(*per_chr[c]):
                    ax.axvspan(s_ / scale, e_ / scale, color=STRIP_COLOURS[cls], lw=0)
                legend_items.append((cls, STRIP_COLOURS[cls]))
        sat.set_ylabel("CenSat / rDNA", rotation=0, ha="right", va="center", fontsize=7)
        sd.set_ylabel("segdups", rotation=0, ha="right", va="center", fontsize=7)

        for ax, clf in zip(axes[n_top:], a.classifiers):
            r = reads[(reads["classifier"] == clf) & (reads["ref"] == c)]
            h = np.histogram(r["mid"], bins=edges)[0]
            for left, n in zip(edges[:-1], h):
                if n > 0:
                    ax.axvspan(left / scale, (left + bin_bp) / scale, ymin=0.1, ymax=0.9,
                               color=cmap(0.15 + 0.85 * norm(n)), lw=0)
            ax.set_ylabel(f"{clf}\n(n={len(r)})", rotation=0, ha="right", va="center", fontsize=8)
            if len(r) == 0:
                ax.text(0.5, 0.55, "no leaked reads", transform=ax.transAxes, ha="center", va="center",
                        fontsize=8, color="#999999")
        for ax in axes:
            ax.set_xlim(0, ln / scale)
            if ax not in axes[:1] or not have_rep:
                ax.set_ylim(0, 1)
            ax.set_yticks([])
            for sp in ("top", "right", "left"):
                ax.spines[sp].set_visible(False)
            if ax in list(axes[:n_top]):
                ax.set_facecolor("#f4f4f4")
        axes[-1].set_xlabel(f"position ({unit})")
        seen = dict(legend_items)
        handles = [plt.Rectangle((0, 0), 1, 1, color=col) for col in seen.values()]
        axes[0].legend(handles, list(seen), loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=6, frameon=False)
        axes[0].set_title(f"{a.title + ': ' if a.title else ''}{rev[c]} leaked reads (pooled: {', '.join(a.samples)})",
                          fontsize=9, loc="left")
        fig.tight_layout(rect=(0, 0, 0.93, 1))
        cax = fig.add_axes([0.945, 0.06, 0.012, 0.26])
        sm = matplotlib.cm.ScalarMappable(norm=mcolors.Normalize(vmin=1, vmax=max(gmax, 2)),
                                          cmap=mcolors.ListedColormap(cmap(np.linspace(0.15, 1, 256))))
        cb = fig.colorbar(sm, cax=cax)
        cb.set_label(f"reads per {bin_bp / scale:g} {unit} bin", fontsize=7)
        cb.ax.tick_params(labelsize=6)
        fig.savefig(os.path.join(a.outdir, f"{rev[c]}.png"), dpi=150)
        pdf.savefig(fig)
        plt.close(fig)
    pdf.close()
    print("wrote", a.outdir)

    # Enrichment: share of leaked-read aligned bases in each annotation vs the share of the genome it covers
    genome_chroms = [c for c in rev if c.startswith("NC_0609") or c == MITO]
    genome = sum(length_of(c) for c in genome_chroms)
    layers = {}
    for src in (repeats, genes, strips):
        layers.update(src)
    cov_genome = {k: sum(int((v[c][1] - v[c][0]).sum()) for c in v if c in set(genome_chroms)) for k, v in layers.items()}
    rows = []
    for clf in a.classifiers:
        r = reads[reads["classifier"] == clf]
        total_bp = int((r["end"] - r["start"] + 1).sum())
        for layer, per_chr in layers.items():
            bp_in = 0
            ge50 = 0
            for c, g in r.groupby("ref"):
                if c not in per_chr:
                    continue
                lo = (g["start"].to_numpy() - 1).astype(np.int64)
                hi = g["end"].to_numpy().astype(np.int64)
                ov = overlap(per_chr[c], lo, hi)
                bp_in += int(ov.sum())
                ge50 += int((ov >= 0.5 * (hi - lo)).sum())
            pg = cov_genome[layer] / genome
            pr = bp_in / total_bp if total_bp else float("nan")
            rows.append({"classifier": clf, "annotation": layer, "n_reads": len(r), "reads_ge50pct_in": ge50,
                         "prop_leaked_bp": pr, "prop_genome": pg, "enrichment": pr / pg if pg and total_bp else float("nan")})
        if "any repeat" in layers and total_bp:
            pr_any = [x for x in rows if x["classifier"] == clf and x["annotation"] == "any repeat"][0]["prop_leaked_bp"]
            pg_any = cov_genome["any repeat"] / genome
            rows.append({"classifier": clf, "annotation": "no repeat (unique)", "n_reads": len(r), "reads_ge50pct_in": "",
                         "prop_leaked_bp": 1 - pr_any, "prop_genome": 1 - pg_any,
                         "enrichment": (1 - pr_any) / (1 - pg_any)})
    pd.DataFrame(rows).to_csv(os.path.join(a.outdir, "enrichment.tsv"), sep="\t", index=False)
    print("wrote enrichment.tsv")


if __name__ == "__main__":
    main()
