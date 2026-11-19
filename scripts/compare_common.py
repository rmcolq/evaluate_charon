#!/usr/bin/env python3
"""Compare classifiers on the reads they all processed, with 95% Wilson intervals.

Charon silently drops reads below its minimum quality / length, which the other tools keep. Counting those reads
for the other tools but not for Charon makes the comparison unfair, so this reports, per sample:

  * all reads each classifier saw
  * "filtered": reads some classifier saw but Charon (or any `--reference`) did not process
  * the host / microbial / unclassified calls restricted to the reads every classifier processed

Reads `<outdir>/<sample>/<sample>_<classifier>_calls.csv.gz`. Writes common_set_comparison.tsv.

  scripts/compare_common.py --outdir out --samples HG02583 HG02922 NA18983 \
      --classifiers charon deacon deacon_sensitive nohuman
"""
import argparse
import math
import os

import duckdb


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - h) / d), min(1.0, (c + h) / d))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="out")
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--classifiers", nargs="+", required=True)
    ap.add_argument("--results", default="results/common_set")
    args = ap.parse_args()
    os.makedirs(args.results, exist_ok=True)

    con = duckdb.connect()
    rows = []
    for sample in args.samples:
        loaded = []
        for clf in args.classifiers:
            path = os.path.join(args.outdir, sample, f"{sample}_{clf}_calls.csv.gz")
            if not os.path.exists(path):
                print(f"skip {path}: missing")
                continue
            con.execute(f"CREATE OR REPLACE TEMP TABLE c_{clf} AS SELECT read_id, CASE WHEN status = 'U' OR classification IS NULL THEN 'unclassified' ELSE classification END AS classification "
                        f"FROM read_csv_auto('{path}')")
            loaded.append(clf)
        if not loaded:
            continue
        union = " UNION ".join(f"SELECT read_id FROM c_{c}" for c in loaded)
        con.execute(f"CREATE OR REPLACE TEMP TABLE all_reads AS {union}")
        common = " INTERSECT ".join(f"SELECT read_id FROM c_{c}" for c in loaded)
        con.execute(f"CREATE OR REPLACE TEMP TABLE common AS {common}")
        n_all = con.execute("SELECT count(*) FROM all_reads").fetchone()[0]
        n_common = con.execute("SELECT count(*) FROM common").fetchone()[0]
        for clf in loaded:
            seen = con.execute(f"SELECT count(*) FROM c_{clf}").fetchone()[0]
            counts = dict(con.execute(
                f"SELECT classification, count(*) FROM c_{clf} JOIN common USING (read_id) GROUP BY 1").fetchall())
            row = {"sample": sample, "classifier": clf, "reads_in_any": n_all, "reads_seen": seen,
                   "not_processed": n_all - seen, "common_reads": n_common}
            for cls in ("human", "microbial", "unclassified"):
                k = counts.get(cls, 0)
                lo, hi = wilson(k, n_common)
                row.update({f"n_{cls}": k, f"prop_{cls}": k / n_common if n_common else float("nan"),
                            f"{cls}_lo": lo, f"{cls}_hi": hi})
            rows.append(row)
    out = os.path.join(args.results, "common_set_comparison.tsv")
    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(out, sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 30):
        print(df[["sample", "classifier", "reads_seen", "not_processed", "common_reads", "n_human", "n_microbial",
                  "n_unclassified", "microbial_lo", "microbial_hi"]].to_string(index=False))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
