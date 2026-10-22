# evaluate_charon
Nextflow pipeline to run charon and collect statistics


## Prebuilt minimap2 index

Every mapping task builds the minimap2 index of the reference from the FASTA, which needs several times
more memory than mapping (about 12 GB for a human reference) and repeats the work for every task. Build the
index once and reuse it:

```
nextflow run main.nf -entry build_mmi_index --refs <reference.fa.gz> --outdir <dir> -profile <profile>
```

This writes `<reference>.map-ont.mmi` (or `.sr.mmi` with `--read_type illumina`) to `--outdir`. Then run the
pipeline with `--refs_mmi <reference>.map-ont.mmi`; it is used instead of `--refs` for all mapping steps. The
index has to be built with the same preset the pipeline maps with, which the entry workflow takes from
`--read_type`. Mapping to the index gives the same alignments as mapping to the FASTA.

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
