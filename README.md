# evaluate_charon
Nextflow pipeline to run charon and collect statistics


## minimap2 index of the reference

Mapping to `--refs` uses a minimap2 index which is built once and kept in `<store_dir>/minimap2/`, so that mapping tasks
load it instead of rebuilding it from the FASTA. The first run which maps (`--evaluate_host` or `--evaluate_microbial`)
builds it, later runs and all methods reuse it, and the path is printed at the start of the run. The size of the
reference is part of the name, so a changed reference is indexed again. `--store_dir` should be an absolute path.

Building needs far more memory than mapping (about 12 GB for a human reference), and the index step runs in the same
container with the same memory limit as other steps. On a machine without that memory, such as a laptop, build the index
natively instead and copy it to the logged path:

```
minimap2 -x map-ont -t 8 -d <store_dir>/minimap2/<name from the log> <reference.fa.gz>
```

The preset has to match the one the pipeline maps with, `map-ont` or `sr` with `--read_type illumina`.

## Leak tiers

Reads called microbial are assigned a `leak_tier` in the per-sample outputs (`*_leaks.parquet`, `*_full.csv`).
Each tier is a subset of the one above:

| Tier | Meaning |
| --- | --- |
| 1 | Called microbial |
| 2 | Tier 1 and maps to the human (T2T) reference, i.e. not unmapped, EBV or another contaminant |
| 3 | Tier 2 and a BLAST human hit with `--blast_min_cov` query coverage and `--blast_min_pident` identity |
| 4 | Tier 3 and the best human bitscore is at least `--blast_bitscore_ratio` times the best non-human bitscore (human-specific, not a shared region) |

Tier 4 needs `--evaluate_microbial` with `--blast_nonhuman` (the default), which adds a BLAST search with human (taxid 9606)
excluded and so needs a local `--blast_db` that includes the taxonomy files.

Real contaminants such as Mycoplasma should be added to the reference in `--refs`, and their accessions listed (whitespace or
comma separated) in a file passed with `--contaminant_accs`, so that they are not scored as human leaks.

Per sample and classifier the pipeline writes `*_calls.parquet` (every read and its classification) and `*_leaks.parquet`
(every read called microbial, with features, mapping, BLAST and tier columns). If `pyarrow` is not in the container, gzipped
CSV is written instead. Run time and peak memory per process are in `execution_trace.txt`.

## Classifiers compared

`--charon`, `--deacon` and `--nohuman` each run one tool with its default settings. `--deacon_sensitive` runs deacon again as a
separate comparator (named `deacon_sensitive` in all outputs) with its high sensitivity settings, an absolute threshold of 1 and
no relative threshold (`-a 1 -r 0`). nohuman has no equivalent mode.

`--charon_host_only` runs charon again with an index of the host only (`--db_host_only`), so reads are called host or not host. Reads
which are not called host are kept, as for deacon and nohuman, and it is named `charon_host_only` in all outputs. It needs a charon
image which supports a host-only index (`--charon_host_only_container`); `--charon_container` sets the image for the normal index.
The host-only index should be built from the same host sequences and k-mer settings as the full index, so that the only difference
between the two is the microbial category.
