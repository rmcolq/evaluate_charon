# evaluate_charon
Nextflow pipeline to run charon and collect statistics

## Paired reads (Illumina)

Add `--paired true` and give paired input, one of:

| Input | Option |
| --- | --- |
| BAM/CRAM of read pairs (e.g. 1000 Genomes high-coverage CRAM) | `--bam_csv` (`sample_id,filepath`), plus `--cram_ref <reference.fa>` for a CRAM |
| FASTQ pairs | `--fastq_pairs_csv` (`sample_id,fastq_1,fastq_2`), or `--fastq` with `--fastq2` and `--unique_id` |

A BAM/CRAM is converted to R1/R2 FASTQ (primary alignments, reads without a mate dropped). Paths can be local or
`s3://`. Each pair travels through the pipeline as one item: charon, deacon and nohuman take both files and call
the pair together, and minimap2 maps pairs with the short-read `sr` preset (the stored minimap2 index is named by
preset, so paired and `--read_type illumina` runs build and keep their own). `--host_eval_reads` counts pairs. In
the evaluation a pair is described by its longer alignment, and `mates_mapped` records how many mates mapped.
Mixed paired and single-end input in one run is not supported, and `--fastq_dir` does not take pairs.

## Build and store minimap2 index of the reference

Mapping to `--refs` uses a minimap2 index which is built once and kept in `<store_dir>/minimap2/`, so that mapping tasks
load it instead of rebuilding it from the FASTA. The first run which maps (`--evaluate_host` or `--evaluate_microbial`)
builds it, later runs and all methods reuse it, and the path is printed at the start of the run. The size of the
reference is part of the name, so a changed reference is indexed again. `--store_dir` should be an absolute path.

Building a human-sized index natively takes about 2 minutes and 3.5 GB of memory (minimap2 2.26, 8 threads). The index step
runs in a container like the other steps, and under emulation or a small memory limit, such as an amd64 container on Apple
Silicon, it can be very slow or be killed. In that case build the index natively and copy it to the logged path:

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
excluded.

BLAST searches the reads which were called microbial and map to the human reference. With `--blast_db` (a directory holding
`nt` and the taxonomy files) the search is local. Without it the search is remote, against NCBI, which is slow and
rate limited, so the pipeline keeps requests small and tolerant of failures:

- reads are sent in chunks of `--blast_chunk_reads` (default 20), one chunk at a time for a remote search;
- each search is tried up to five times with a growing wait;
- `--blast_max_reads` (default: all) searches a fixed random sample of that many reads per sample and method;
- results are kept per read in `<store_dir>/blast_cache/<sample>` (`--blast_cache false` keeps them in the task instead), so a
  read which leaked in several methods, or in an earlier run, is only searched once. `--store_dir` should be an absolute path
  on a filesystem every task can see;
- a read whose search still fails is written as `failed` and gets `blast_failed = true` in the leak tables. It is counted in
  `num_microbial_map_host_blast_failed_<classifier>` and is neither verified nor shown to be non-human. The run does not stop and
  the summary is still written.

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
