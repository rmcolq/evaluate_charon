#!/usr/bin/env python3
"""Compare Charon threshold combinations across datasets, from the combination_grid.tsv each
charon_param_sweep.py run writes. Writes one wide table and a plot of human reads called microbial
against microbial reads kept, one panel per dataset, with the run's settings and any --highlight settings marked.

    scripts/compare_param_sweeps.py --sweep sim_ont=sim/param_sweep art_ont=art/param_sweep \\
        --highlight 7,0.07,0.04,0,0 --outdir comparison
A --highlight is confidence,min_proportion_diff,host_unique_prop_lo_threshold,min_probability_diff,min_compression.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

KEYS = ["confidence_threshold", "min_proportion_difference", "host_unique_prop_lo_threshold",
        "min_prob_difference", "min_compression"]
RUN = (7, 0.04, 0.05, 0, 0)
MARKS = ["#eb6834", "#1baf7a", "#e87ba4", "#eda100", "#008300"]
INK, MUTED, GRID, POINTS = "#2b2b29", "#6b6a64", "#e4e3de", "#9ec5f4"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", nargs="+", required=True, help="name=dir pairs, each dir holding combination_grid.tsv")
    ap.add_argument("--highlight", nargs="*", default=[], help="settings to mark, as comma-separated values of the five keys")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                         "grid.color": GRID, "axes.axisbelow": True, "legend.frameon": False, "figure.dpi": 150})

    sweeps = dict(s.split("=", 1) for s in a.sweep)
    highlights = [tuple(float(v) for v in h.split(",")) for h in a.highlight]
    wide = None
    fig, axes = plt.subplots(1, len(sweeps), figsize=(5 * len(sweeps), 4.5), constrained_layout=True, squeeze=False)
    for ax, (name, d) in zip(axes[0], sweeps.items()):
        g = pd.read_csv(os.path.join(d, "combination_grid.tsv"), sep="\t")
        cols = {"human_microbial": "human_called_microbial", "microbial_microbial_pct": "microbial_kept_pct",
                "microbial_microbial_bases_pct": "microbial_kept_bases_pct", "human_human_pct": "human_removed_pct",
                "microbial_human": "microbial_called_human"}
        g = g[KEYS + list(cols)].rename(columns=cols)
        ax.scatter(g.microbial_kept_pct, g.human_called_microbial, s=8, color=POINTS, lw=0, label="all combinations")
        for setting, color, label in [(RUN, INK, "run settings")] + \
                [(h, MARKS[i % len(MARKS)], ",".join(f"{v:g}" for v in h)) for i, h in enumerate(highlights)]:
            row = g[(g[KEYS].sub(setting).abs() < 1e-9).all(axis=1)]
            if len(row):
                ax.scatter(row.microbial_kept_pct, row.human_called_microbial, s=60, facecolor="none",
                           edgecolor=color, lw=1.8, label=label)
        ax.set_title(name, loc="left")
        ax.set_xlabel("Microbial reads called microbial (%)")
        ax.set_ylabel("Human reads called microbial")
        ax.legend(fontsize=7)
        g = g.rename(columns={c: f"{name}:{c}" for c in cols.values()})
        wide = g if wide is None else wide.merge(g, on=KEYS)
    fig.savefig(os.path.join(a.outdir, "combinations_by_dataset.png"))
    wide.to_csv(os.path.join(a.outdir, "combinations_by_dataset.tsv"), sep="\t", index=False)
    print(wide[(wide[KEYS].sub(RUN).abs() < 1e-9).all(axis=1)].T.to_string())


if __name__ == "__main__":
    main()
