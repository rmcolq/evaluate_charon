#!/usr/bin/env python3
"""Pool the BLAST queries of every sample and method into one set of chunks.

Each search reads the whole BLAST database, so the time is set by the number of searches, not the number of reads.
The reads of all samples and methods are therefore searched together: a read is searched once per sample, however
many methods let it through, and reads already in the cache are not searched again.

Input: a manifest (sample, label) with one line per query file, in the order of the files in --queries (named query,
or query1, query2, ... as Nextflow stages them); a sample and label may have several files. Files may be gzipped.
The label is, for example, the method whose reads these are.
Output:
- chunk_NNNN.fasta: the reads to search, named q<N> (read ids of different samples may be the same);
- map.tsv: chunk, q<N>, sample, read id;
- selected.tsv: sample, method, read id, for every read whose result the sample and method need (after any cap);
- groups.tsv: sample, method, for every sample and method, including those with no reads;
- published/<sample>/<sample>_<method>_blast_queries.fasta.gz: the reads of each sample and method, before any cap.
"""
import argparse
import gzip
import os
import random
import re


def safe(read_id):
    """The cache file name of a read (the same rule as before pooling, so existing caches still match)."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", read_id)


def open_text(path):
    with open(path, "rb") as handle:
        gzipped = handle.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gzipped else open(path)


def staged_files(directory):
    """The staged query files in their original order: query, or query1, query2, ..."""
    names = os.listdir(directory)
    return [os.path.join(directory, name) for name in sorted(names, key=lambda n: int(n[len("query"):] or 0))]


def read_fasta(path):
    """(id, sequence) for each record; the id is the header up to the first space or tab; wrapped FASTA is fine."""
    read_id, seq = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if read_id is not None:
                    yield read_id, "".join(seq)
                read_id, seq = re.split(r"[ \t]", line[1:], maxsplit=1)[0], []
            elif line:
                seq.append(line)
    if read_id is not None:
        yield read_id, "".join(seq)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True, help="TSV: sample, label; one line per file in --queries")
    parser.add_argument("--queries", required=True, help="directory of the staged query files")
    parser.add_argument("--cache", default="", help="cache directory (<cache>/<sample>/<read>.done); empty for none")
    parser.add_argument("--max_reads", type=int, default=0, help="search at most this many reads per sample and method (0: all)")
    parser.add_argument("--chunk_reads", type=int, required=True, help="reads per chunk")
    args = parser.parse_args()

    # unique reads of each sample and method, in the order first seen
    reads = {}  # (sample, method) -> {read id: sequence}
    with open(args.manifest) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    files = staged_files(args.queries)
    if len(files) != len(rows):
        raise SystemExit(f"{len(rows)} manifest lines but {len(files)} query files")
    for (sample, method), fasta in zip(rows, files):
        group = reads.setdefault((sample, method), {})
        for read_id, seq in read_fasta(fasta):
            group.setdefault(read_id, seq)

    os.makedirs("published", exist_ok=True)
    selected = {}  # (sample, method) -> [read id]
    for (sample, method), group in reads.items():
        os.makedirs(f"published/{sample}", exist_ok=True)
        with gzip.open(f"published/{sample}/{sample}_{method}_blast_queries.fasta.gz", "wt") as out:
            for read_id, seq in group.items():
                out.write(f">{read_id}\n{seq}\n")
        ids = list(group)
        if args.max_reads and len(ids) > args.max_reads:
            # a fixed random sample: the same reads every run, so cached results are reused
            keep = set(random.Random(f"{sample}\t{method}").sample(ids, args.max_reads))
            ids = [read_id for read_id in ids if read_id in keep]
        selected[(sample, method)] = ids

    with open("groups.tsv", "w") as out:
        for sample, method in reads:
            out.write(f"{sample}\t{method}\n")
    with open("selected.tsv", "w") as out:
        for (sample, method), ids in selected.items():
            for read_id in ids:
                out.write(f"{sample}\t{method}\t{read_id}\n")

    # each read of a sample once, unless the cache already has its result
    todo = {}  # (sample, read id) -> sequence
    for (sample, method), ids in selected.items():
        for read_id in ids:
            key = (sample, read_id)
            if key in todo:
                continue
            if args.cache and os.path.exists(os.path.join(args.cache, sample, safe(read_id) + ".done")):
                continue
            todo[key] = reads[(sample, method)][read_id]

    chunk_reads = max(1, args.chunk_reads)
    with open("map.tsv", "w") as mapping:
        chunk = None
        for n, ((sample, read_id), seq) in enumerate(todo.items()):
            if n % chunk_reads == 0:
                if chunk:
                    chunk.close()
                name = f"chunk_{n // chunk_reads:04d}"
                chunk = open(f"{name}.fasta", "w")
            chunk.write(f">q{n}\n{seq}\n")
            mapping.write(f"{name}\tq{n}\t{sample}\t{read_id}\n")
        if chunk:
            chunk.close()
    print(f"{len(todo)} reads to search in {-(-len(todo) // chunk_reads)} chunks "
          f"({sum(len(ids) for ids in selected.values())} reads needed by {len(selected)} samples and methods)")


if __name__ == "__main__":
    main()
