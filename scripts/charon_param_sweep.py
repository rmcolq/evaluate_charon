#!/usr/bin/env python3
"""Re-call Charon dehost results under different thresholds and score them against simulated truth.

Every threshold applied in ReadEntry::call_host (include/read_entry.hpp) uses values printed in the
.out/.tsv report (confidence, length, mean quality, compression, and per-category unique proportion
and probability), so for reads that were classified (status C or U) the call under any other threshold
is exact. Reads with status F were filtered before hashing, so lowering min_length or min_quality below
the values used in the run cannot be simulated. Parameters that change the fitted model
(lo_hi_threshold, num_reads_to_fit, dist) need a rerun.

Truth comes from simulated read headers (`--headers`: Badread, or an accession in the read id) or from a
read_id<TAB>label table (`--truth_tsv`, e.g. read2truth.ont.tsv for the Hall et al. real-read mixture).
Pass the run's log (`--log`) so the thresholds it used are the starting point.
"""
import argparse
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HUMAN_PREFIXES = ("CM0356", "CM0357", "JAHRJT")  # KOREF chromosomes and unplaced contigs
MTB = {"NZ_OW052570.1"}

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"
INK, MUTED, GRID = "#2b2b29", "#6b6a64", "#e4e3de"

DEFAULTS = dict(confidence_threshold=7, min_proportion_difference=0.04, host_unique_prop_lo_threshold=0.05,
                min_prob_difference=0.0, confidence_probability_threshold=0.0,
                min_length=140, min_quality=15.0, min_compression=0.0)


def style():
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                         "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
                         "lines.linewidth": 2, "legend.frameon": False, "figure.dpi": 150})


def classify_accessions(segs, human_prefixes, mtb):
    """truth and source for the source accessions of one read (more than one for a Badread chimera)"""
    human = [a.split("|")[-1].startswith(human_prefixes) for a in segs]
    if all(human):
        return "human", "human"
    if any(human):
        return "chimera", "chimera"
    accs = {a.split("|")[-1] for a in segs}
    if accs <= mtb:
        return "microbial", "MTB"
    if not any("|" in a for a in segs):
        return "microbial", "NTM"
    return "microbial", "other microbial"


def load_truth_headers(path, human_prefixes=HUMAN_PREFIXES, mtb=MTB):
    """Truth from simulated read headers, one per line.

    Badread:            @<id> <acc>,<strand>,<start>-<end> [chimera <acc>,...] length=...
    accession in the id: @<acc>-<n>/1 (simulated Illumina pairs)
    Lines that are not headers (a FASTQ quality line starting with @) find no report row and drop out of the join.
    """
    rows = []
    with open(path) as f:
        for line in f:
            parts = line[1:].split()
            if not parts:
                continue
            segs = [p.split(",")[0] for p in parts[1:] if "strand" in p]
            if not segs:
                segs = [parts[0].rsplit("/", 1)[0].rsplit("-", 1)[0]]
            truth, source = classify_accessions(segs, human_prefixes, mtb)
            rows.append((parts[0], truth, source, segs[0].split("|")[-1]))
    return pd.DataFrame(rows, columns=["id", "truth", "source", "accession"])


def load_truth_tsv(path, human_labels=("human",)):
    """Truth from a two-column read_id<TAB>label table (e.g. read2truth.ont.tsv); labels other than human are microbial"""
    t = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["id", "source"], dtype=str)
    t["truth"] = np.where(t.source.isin(human_labels), "human", "microbial")
    return t


DETAIL_RE = r"(?P<cat>[^: ]+):(?P<count>[^: ]+):(?P<prop>[^: ]+):(?P<uprop>[^: ]+):(?P<prob>[^: ]+)"


def load_report(path, host="human"):
    r = pd.read_csv(path, sep="\t", header=None, usecols=range(9), quoting=3, dtype={"details": str, "call": str},
                    names=["status", "id", "call", "length", "num_hashes", "quality", "confidence", "compression", "details"])
    r["call"] = r["call"].fillna("")
    r.loc[r.status == "U", "call"] = "unclassified"
    r.loc[r.status == "F", "call"] = "filtered"
    # details: "cat:count:prop:uprop:prob cat:count:prop:uprop:prob" (two categories for a host + microbial index)
    det = r.details.where(r.status != "F", "").fillna("").str.extractall(DETAIL_RE)
    det = det.apply(lambda c: c if c.name == "cat" else pd.to_numeric(c, errors="coerce"))
    is_host = det.cat == host
    for prefix, sel in (("host", is_host), ("other", ~is_host)):
        part = det[sel].droplevel(1)
        for col in ("count", "prop", "uprop", "prob"):
            r[f"{prefix}_{col}"] = part[col].reindex(r.index)
    r["confidence"] = pd.to_numeric(r.confidence, errors="coerce")
    r["compression"] = pd.to_numeric(r.compression, errors="coerce")
    return r


LOG_KEYS = {"min_length": "min_length", "min_quality": "min_quality", "min_compression": "min_compression",
            "confidence_threshold": "confidence_threshold",
            "confidence_probability_threshold": "confidence_probability_threshold",
            "host_unique_prop_lo_threshold": "host_unique_prop_lo_threshold",
            "min_proportion_difference": "min_proportion_difference", "min_prob_difference": "min_prob_difference"}


def thresholds_from_log(path):
    """The thresholds a dehost run used, from the 'Dehost Arguments' block of its log"""
    found = {}
    with open(path, encoding="utf-8-sig") as f:
        for line in f:
            key, _, value = line.strip().partition(":")
            if key in LOG_KEYS and value.strip():
                found[LOG_KEYS[key]] = float(value)
    return found


def recall(df, **overrides):
    """Vectorised port of ReadEntry::call_host. Returns array of human / microbial / unclassified / filtered."""
    p = dict(DEFAULTS, **overrides)
    a = df.attrs["arrays"]
    filtered = a["F"].copy()
    if p["min_length"] > DEFAULTS["min_length"] or p["min_quality"] > DEFAULTS["min_quality"]:
        # raising the pre-filter thresholds: these reads would now be F
        filtered |= (a["length"] < p["min_length"]) | (a["quality"] < p["min_quality"])
    hu, ou, hp, op, conf = a["hu"], a["ou"], a["hp"], a["op"], a["conf"]
    ok = (conf >= p["confidence_threshold"]) & (a["quality"] >= p["min_quality"]) \
        & (a["length"] >= p["min_length"]) & (a["compression"] >= p["min_compression"])
    host = ok & (hu > ou) & (hu - ou > p["min_proportion_difference"]) & (hp > op) \
        & (hp - op > p["min_prob_difference"]) & (np.maximum(hp * conf, conf) >= p["confidence_probability_threshold"])
    micro = ok & ~host & (hu < p["host_unique_prop_lo_threshold"]) & (hu < ou) \
        & (ou - hu > p["min_proportion_difference"]) & (hp < op) & (op - hp > p["min_prob_difference"]) \
        & (np.maximum(op * conf, conf) >= p["confidence_probability_threshold"])
    out = np.full(len(df), 2, dtype=np.int8)  # codes index CALLS
    out[host] = 0
    out[micro] = 1
    out[filtered] = 3
    return out


CALLS = np.array(["human", "microbial", "unclassified", "filtered"], dtype=object)


def prepare(df):
    with np.errstate(invalid="ignore"):
        df.attrs["arrays"] = dict(
            F=(df.status == "F").to_numpy(), length=df.length.to_numpy(float), quality=df.quality.to_numpy(float),
            compression=df.compression.to_numpy(float), conf=df.confidence.to_numpy(float),
            hu=df.host_uprop.to_numpy(float), ou=df.other_uprop.to_numpy(float),
            hp=df.host_prob.to_numpy(float), op=df.other_prob.to_numpy(float),
            truth={t: (df.truth == t).to_numpy() for t in ("human", "microbial", "chimera")})
    return df


def score(df, codes):
    a = df.attrs["arrays"]
    res = {}
    for truth, m in a["truth"].items():
        n = m.sum()
        reads = np.bincount(codes[m], minlength=4)
        bases = np.bincount(codes[m], weights=a["length"][m], minlength=4)
        total_bases = bases.sum()
        for k, c in enumerate(CALLS):
            res[f"{truth}_{c}"] = int(reads[k])
            res[f"{truth}_{c}_pct"] = 100 * reads[k] / n if n else np.nan
            res[f"{truth}_{c}_bases_pct"] = 100 * bases[k] / total_bases if n else np.nan
    res["human_kept"] = res["human_microbial"] + res["human_unclassified"] + res["human_filtered"]
    return res


SWEEPS = {
    "confidence_threshold": np.r_[0:41, 50, 60, 80, 100, 150, 200, 255],
    "min_proportion_difference": np.round(np.r_[0:0.3001:0.005], 4),
    "host_unique_prop_lo_threshold": np.round(np.r_[0:0.5001:0.01, 0.6, 0.8, 1.0], 4),
    "min_prob_difference": np.r_[0, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999],
    "min_length": np.r_[80, 100, 120, 140, 160, 200, 250, 300, 500, 750, 1000, 1500, 2000, 3000, 5000],
    "min_quality": np.round(np.r_[0:25.01:0.5], 2),
    "min_compression": np.round(np.r_[0:0.4001:0.01], 3),
}
LABELS = {
    "confidence_threshold": "--confidence", "min_proportion_difference": "--min_proportion_diff",
    "host_unique_prop_lo_threshold": "--host_unique_prop_lo_threshold", "min_prob_difference": "--min_probability_diff",
    "min_length": "--min_length", "min_quality": "--min_quality", "min_compression": "--min_compression",
}


XLIM = {"confidence_threshold": (0, 60), "host_unique_prop_lo_threshold": (0, 0.2),
        "min_compression": (0, 0.32), "min_proportion_difference": (0, 0.2), "min_quality": (15, 25)}


def plot_sweeps(sweep, out):
    params = list(SWEEPS)
    fig, axes = plt.subplots(len(params), 3, figsize=(12, 2.3 * len(params)), constrained_layout=True)
    for i, prm in enumerate(params):
        s = sweep[sweep.parameter == prm]
        if prm == "min_prob_difference":  # spans 0..0.999: plot on an index axis with value labels
            x = np.arange(len(s))
            default_x = 0
        else:
            x = s.value.astype(float).to_numpy()
            default_x = DEFAULTS[prm]
        human_not_removed = 100 * (1 - s.human_human / (s.human_human + s.human_microbial + s.human_unclassified + s.human_filtered))
        panels = [(s.human_microbial, ORANGE, "human called microbial"),
                  (human_not_removed, YELLOW, "human not called human"),
                  (s.microbial_microbial_pct, BLUE, "microbial called microbial")]
        for a, (y, color, lab) in zip(axes[i], panels):
            a.plot(x, y, color=color, marker="o", ms=2.5)
            a.axvline(default_x, color=MUTED, ls="--", lw=1)
            a.ticklabel_format(axis="y", useOffset=False, style="plain")
            if prm in XLIM:
                lo, hi = XLIM[prm]
                a.set_xlim(lo, hi)
                vis = y[(x >= lo) & (x <= hi)]
                pad = 0.05 * (vis.max() - vis.min() + 1e-9)
                a.set_ylim(vis.min() - pad, vis.max() + pad)
            if prm == "min_prob_difference":
                a.set_xticks(x[::2], [f"{v:g}" for v in s.value.to_numpy()[::2]], rotation=45, fontsize=7)
        axes[i, 0].set_ylabel(LABELS[prm], fontsize=9, fontweight="bold")
        if i == 0:
            axes[i, 0].set_title("Human reads called microbial (count)", loc="left")
            axes[i, 1].set_title("Human reads not called human (%)", loc="left")
            axes[i, 2].set_title("Microbial reads called microbial (%)", loc="left")
    fig.suptitle("Charon dehost on the simulated metagenome, one parameter changed at a time (dashed = value in the run).\n"
                 "No microbial read is called human at any setting shown.", x=0.01, ha="left", fontsize=10)
    fig.savefig(out)
    plt.close(fig)


def plot_tradeoff(sweep, out):
    colors = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, MUTED]
    fig, a = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for prm, color in zip(SWEEPS, colors):
        s = sweep[sweep.parameter == prm].sort_values("microbial_microbial_pct", ascending=False)
        lost = 100 - s.microbial_microbial_pct
        keep = lost <= 3
        a.plot(lost[keep], s.human_microbial[keep], color=color, marker="o", ms=3, label=LABELS[prm], alpha=0.9)
    b = sweep[(sweep.parameter == "confidence_threshold") & (sweep.value == DEFAULTS["confidence_threshold"])].iloc[0]
    a.plot(100 - b.microbial_microbial_pct, b.human_microbial, "o", ms=10, mfc="none", mec=INK, mew=1.5)
    a.annotate("run settings", (100 - b.microbial_microbial_pct, b.human_microbial), xytext=(8, 6),
               textcoords="offset points", fontsize=8, color=INK)
    a.set_xlabel("Microbial reads not called microbial (%)")
    a.set_ylabel("Human reads called microbial (count)")
    a.set_title("Trade-off for each parameter moved on its own (points within 3% microbial loss)", loc="left")
    a.legend(fontsize=8)
    fig.savefig(out)
    plt.close(fig)


def plot_distributions(df, out):
    d = df[df.status != "F"].copy()
    d["uprop_diff"] = d.host_uprop - d.other_uprop
    d["prob_diff"] = d.host_prob - d.other_prob
    feats = [("confidence", "confidence (capped at 255)", DEFAULTS["confidence_threshold"], np.r_[0:256:3]),
             ("uprop_diff", "host − microbial unique proportion", [-DEFAULTS["min_proportion_difference"], DEFAULTS["min_proportion_difference"]], np.r_[-1:1.001:0.02]),
             ("host_uprop", "host unique proportion", DEFAULTS["host_unique_prop_lo_threshold"], np.r_[0:1.001:0.01]),
             ("prob_diff", "host − microbial probability", None, np.r_[-1:1.001:0.02]),
             ("length", "read length (log10)", np.log10(DEFAULTS["min_length"]), np.r_[2:5.5:0.05]),
             ("quality", "mean quality", DEFAULTS["min_quality"], np.r_[0:40:0.5]),
             ("compression", "compression", None, np.r_[0:0.6:0.005])]
    fig, axes = plt.subplots(len(feats), 1, figsize=(9, 2.1 * len(feats)), constrained_layout=True)
    for a, (col, lab, thr, bins) in zip(axes, feats):
        for truth, color in (("human", ORANGE), ("microbial", BLUE), ("chimera", MAGENTA)):
            v = d.loc[d.truth == truth, col]
            if v.empty:
                continue
            if col == "length":
                v = np.log10(v)
            a.hist(v, bins=bins, histtype="step", lw=1.5, color=color, label=f"truth {truth}", density=False)
        a.set_yscale("log")
        a.set_xlabel(lab)
        for t in np.atleast_1d(thr) if thr is not None else []:
            a.axvline(t, color=MUTED, ls="--", lw=1)
    axes[0].legend(fontsize=8)
    axes[0].set_title("Report values by true origin (reads not filtered; dashed = threshold in run)", loc="left")
    fig.savefig(out)
    plt.close(fig)


def plot_grid(df, out):
    confs = np.r_[0, 2, 4, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50]
    mpds = np.round(np.r_[0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1, 0.12, 0.15, 0.2], 3)
    leak = np.zeros((len(confs), len(mpds)))
    ret = np.zeros_like(leak)
    for i, c in enumerate(confs):
        for j, m in enumerate(mpds):
            s = score(df, recall(df, confidence_threshold=c, min_proportion_difference=m))
            leak[i, j] = s["human_microbial"]
            ret[i, j] = s["microbial_microbial_pct"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for a, z, title, cmap, fmt in ((axes[0], leak, "Human reads called microbial", "Oranges", "{:.0f}"),
                                   (axes[1], ret, "Microbial reads called microbial (%)", "Blues", "{:.2f}")):
        im = a.imshow(z, cmap=cmap, aspect="auto", origin="lower")
        a.set_xticks(range(len(mpds)), mpds, rotation=45)
        a.set_yticks(range(len(confs)), confs)
        a.set_xlabel("--min_proportion_diff")
        a.set_ylabel("--confidence")
        a.set_title(title, loc="left")
        a.grid(False)
        hi = z.max()
        for i in range(z.shape[0]):
            for j in range(z.shape[1]):
                a.text(j, i, fmt.format(z[i, j]), ha="center", va="center", fontsize=6.5,
                       color="white" if (z[i, j] - z.min()) > 0.6 * (hi - z.min() + 1e-9) else INK)
        a.add_patch(plt.Rectangle((list(mpds).index(0.04) - .5, list(confs).index(7) - .5), 1, 1, fill=False, ec=INK, lw=1.5))
    fig.savefig(out)
    plt.close(fig)
    return pd.DataFrame([(c, m, leak[i, j], ret[i, j]) for i, c in enumerate(confs) for j, m in enumerate(mpds)],
                        columns=["confidence_threshold", "min_proportion_difference", "human_called_microbial", "microbial_retained_pct"])


def plot_sources(df, out):
    order = ["human", "chimera", "MTB", "NTM", "other microbial"]
    calls = ["human", "microbial", "unclassified", "filtered"]
    colors = [ORANGE, BLUE, YELLOW, MUTED]
    tab = pd.crosstab(df.source, df.call).reindex(index=order, columns=calls).fillna(0)
    pct = 100 * tab.div(tab.sum(axis=1), axis=0)
    fig, a = plt.subplots(figsize=(9, 3.2), constrained_layout=True)
    left = np.zeros(len(order))
    for c, col in zip(calls, colors):
        a.barh(order, pct[c], left=left, color=col, label=f"called {c}", edgecolor="white", linewidth=1)
        left += pct[c]
    for k, src in enumerate(order):
        a.text(101, k, f"n={int(tab.loc[src].sum()):,}", va="center", fontsize=8, color=MUTED)
    a.set_xlim(0, 112)
    a.invert_yaxis()
    a.set_xlabel("% of reads")
    a.legend(ncol=4, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=8)
    fig.savefig(out)
    plt.close(fig)
    return tab


COMBO = dict(confidence_threshold=[5, 7, 8, 10],
             min_proportion_difference=[0.04, 0.05, 0.06, 0.07, 0.08, 0.1],
             host_unique_prop_lo_threshold=[0.03, 0.04, 0.05],
             min_prob_difference=[0, 0.2],
             min_compression=[0, 0.25, 0.28, 0.3])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    truth = ap.add_mutually_exclusive_group(required=True)
    truth.add_argument("--headers", help="read headers from the simulated FASTQ (one per line); truth is the source accession")
    truth.add_argument("--truth_tsv", help="read_id<TAB>label table, e.g. read2truth.ont.tsv")
    ap.add_argument("--report", required=True, help="charon dehost report (.tsv/.out)")
    ap.add_argument("--log", help="charon dehost log; the thresholds the run used are read from it (default: Charon defaults)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--human_prefixes", nargs="+", default=list(HUMAN_PREFIXES),
                    help="with --headers: accession prefixes of human source contigs (default: KOREF, as in Hall et al. Zenodo 8339789)")
    ap.add_argument("--mtb", nargs="+", default=sorted(MTB), help="with --headers: accessions reported separately as MTB")
    ap.add_argument("--human_labels", nargs="+", default=["human"], help="with --truth_tsv: labels that are human")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    style()
    if a.log:
        DEFAULTS.update(thresholds_from_log(a.log))
    print("run thresholds:", DEFAULTS)

    t = load_truth_headers(a.headers, tuple(a.human_prefixes), set(a.mtb)) if a.headers \
        else load_truth_tsv(a.truth_tsv, a.human_labels)
    df = prepare(t.merge(load_report(a.report), on="id", how="inner"))
    assert len(df), "no read ids in common between truth and report"

    # check the port reproduces the run
    recalled = CALLS[recall(df)]
    mism = (recalled != df.call.to_numpy())
    print(f"{len(df):,} reads; port reproduces {100 * (1 - mism.mean()):.4f}% of calls ({mism.sum()} differ)")
    if mism.any():
        df.loc[mism].assign(recalled=recalled[mism]).to_csv(os.path.join(a.outdir, "recall_mismatches.tsv"), sep="\t", index=False)

    base = pd.crosstab([df.truth, df.source], df.call, margins=True)
    base.to_csv(os.path.join(a.outdir, "baseline_truth_vs_call.tsv"), sep="\t")
    print(base)

    rows = []
    for prm, vals in SWEEPS.items():
        if prm in ("min_length", "min_quality"):  # lowering these needs a rerun (status F reads were never hashed)
            vals = vals[vals >= DEFAULTS[prm]]
        for v in vals:
            rows.append(dict(parameter=prm, value=v, **score(df, recall(df, **{prm: v}))))
    sweep = pd.DataFrame(rows)
    sweep.to_csv(os.path.join(a.outdir, "one_at_a_time_sweep.tsv"), sep="\t", index=False)

    rows = []
    for combo in itertools.product(*COMBO.values()):
        settings = dict(zip(COMBO, combo))
        rows.append(dict(**settings, **score(df, recall(df, **settings))))
    pd.DataFrame(rows).to_csv(os.path.join(a.outdir, "combination_grid.tsv"), sep="\t", index=False)

    plot_sweeps(sweep, os.path.join(a.outdir, "sweep_one_at_a_time.png"))
    plot_tradeoff(sweep, os.path.join(a.outdir, "tradeoff.png"))
    plot_distributions(df, os.path.join(a.outdir, "feature_distributions.png"))
    plot_grid(df, os.path.join(a.outdir, "grid_confidence_x_proportion_diff.png")) \
        .to_csv(os.path.join(a.outdir, "grid_confidence_x_proportion_diff.tsv"), sep="\t", index=False)

    # human reads by source accession, so leaks concentrated in a few contigs (or a contig set) show up
    if "accession" in df:
        h = df[df.truth == "human"]
        group = h.accession.str.extract(r"^([A-Z_]+)", expand=False)
        pd.crosstab(group, h.call).to_csv(os.path.join(a.outdir, "human_calls_by_accession_prefix.tsv"), sep="\t")
        pd.crosstab(h.accession, h.call).sort_values("microbial", ascending=False) \
            .to_csv(os.path.join(a.outdir, "human_calls_by_accession.tsv"), sep="\t")
    err = df[((df.truth == "human") & (df.call != "human")) | ((df.truth == "microbial") & (df.call != "microbial"))
             | (df.truth == "chimera")]
    err.drop(columns=["details"]).to_csv(os.path.join(a.outdir, "reads_not_called_as_truth.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
