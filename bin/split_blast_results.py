#!/usr/bin/env python3
"""Give each sample and method its BLAST results from the pooled search (see pool_blast_queries.py).

Each result line is "<search>\t<BLAST line>", <search> being "all" or "nonhuman"; the query id (q<N>) is replaced by
the read id. A read whose search failed, or whose chunk has no result file (its task failed), is written as
"failed\t<read id>". New results are added to the cache, if there is one; reads which were not searched this time
take their results from it.

Output: out/<sample>/blast/<sample>.<method>.results_blastn.txt for every sample and method in groups.tsv.
"""
import argparse
import os
import re
from collections import defaultdict


def safe(read_id):
    return re.sub(r"[^A-Za-z0-9._-]", "_", read_id)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--map", required=True)
    parser.add_argument("--selected", required=True)
    parser.add_argument("--groups", required=True)
    parser.add_argument("--results_dir", required=True, help="directory of results_<chunk>.txt files")
    parser.add_argument("--cache", default="", help="cache directory; empty for none")
    args = parser.parse_args()

    query = {}  # q<N> -> (sample, read id)
    chunk_queries = defaultdict(list)
    with open(args.map) as handle:
        for line in handle:
            chunk, qid, sample, read_id = line.rstrip("\n").split("\t")
            query[qid] = (sample, read_id)
            chunk_queries[chunk].append(qid)

    searched = set()  # (sample, read id) with a finished search
    failed = set()
    hits = defaultdict(list)  # (sample, read id) -> result lines
    for chunk, qids in chunk_queries.items():
        path = os.path.join(args.results_dir, f"results_{chunk}.txt")
        if not os.path.exists(path):
            failed.update(query[qid] for qid in qids)
            continue
        chunk_failed = set()
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "failed":
                    chunk_failed.add(fields[1])
                    continue
                sample, read_id = query[fields[1]]
                fields[1] = read_id
                hits[(sample, read_id)].append("\t".join(fields))
        for qid in qids:
            (failed if qid in chunk_failed else searched).add(query[qid])

    if args.cache:
        for sample, read_id in searched:
            base = os.path.join(args.cache, sample, safe(read_id))
            os.makedirs(os.path.dirname(base), exist_ok=True)
            if hits[(sample, read_id)]:
                with open(base + ".txt", "w") as out:
                    out.write("".join(line + "\n" for line in hits[(sample, read_id)]))
            elif os.path.exists(base + ".txt"):
                os.remove(base + ".txt")
            open(base + ".done", "w").close()

    selected = defaultdict(list)
    with open(args.selected) as handle:
        for line in handle:
            sample, method, read_id = line.rstrip("\n").split("\t")
            selected[(sample, method)].append(read_id)

    with open(args.groups) as handle:
        groups = [tuple(line.rstrip("\n").split("\t")) for line in handle if line.strip()]
    for sample, method in groups:
        os.makedirs(f"out/{sample}/blast", exist_ok=True)
        with open(f"out/{sample}/blast/{sample}.{method}.results_blastn.txt", "w") as out:
            for read_id in selected[(sample, method)]:
                key = (sample, read_id)
                if key in searched:
                    out.write("".join(line + "\n" for line in hits[key]))
                    continue
                if key not in failed and args.cache:
                    base = os.path.join(args.cache, sample, safe(read_id))
                    if os.path.exists(base + ".done"):
                        if os.path.exists(base + ".txt"):
                            with open(base + ".txt") as cached:
                                out.write(cached.read())
                        continue
                out.write(f"failed\t{read_id}\n")


if __name__ == "__main__":
    main()
