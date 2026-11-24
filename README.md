# evaluate_charon
Nextflow pipeline to run charon and collect statistics

## Environment

`environment.yml` has Nextflow (with Java) and the Python packages for the analysis scripts in `scripts/`:

```
conda env create -f environment.yml
conda activate evaluate_charon
```

The pipeline's own processes run in containers, so Docker (`-profile docker`) is needed as well. The environment has not
been solved on every platform: if a package fails to resolve, relax its version in the file.

## Paired reads (Illumina)

Give paired input, one of:

| Input | Option |
| --- | --- |
| BAM/CRAM of read pairs (e.g. 1000 Genomes high-coverage CRAM) | `--bam_csv` (`sample_id,filepath`) with `--paired true`, plus `--cram_ref <reference.fa>` for a CRAM |
| FASTQ pairs | `--fastq_pairs_csv` (`sample_id,fastq_1,fastq_2`), or `--fastq` with `--fastq2` (and `--unique_id`; default the first file's name) |

FASTQ pairs are recognised as pairs without `--paired`; a BAM/CRAM needs it, because the file does not say whether it
holds pairs.

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

Tiers 2-4 need the microbial reads mapped, which `--evaluate_microbial` switches on. Without it only tier 1 is
measured: the summary leaves the `*_microbial_unmapped/map_ebv/map_host/map_host_verified` and tier 2-4 columns empty
(they are not zero), and `microbial_mapped` is `False` in the per-read tables. Likewise, without `--evaluate_host` the
`*_host_sampled/unmapped/map_ebv/map_host` columns are empty.

Tier 4 needs `--evaluate_microbial` with `--blast_nonhuman` (the default), which adds a BLAST search with human (taxid 9606)
excluded.

BLAST searches the reads which were called microbial and map to the human reference. With `--blast_db` (a directory holding
`nt` and the taxonomy files) the search is local. Without it the search is remote, against NCBI, which is slow and
rate limited, so the pipeline keeps requests small and tolerant of failures:

- reads are sent in chunks of `--blast_chunk_reads` (default 20), one chunk at a time for a remote search;
- each search is tried up to five times with a growing wait. A remote search is abandoned after `--blast_remote_timeout` (default `20m`)
  times the attempt number (20m, 40m, 60m, 80m, 100m), so a slow search gets longer on later attempts, and a search which
  finished is not repeated;
- `--blast_max_reads` (default: all) searches a fixed random sample of that many reads per sample and method;
- results are kept per read in `<store_dir>/blast_cache/<sample>` (`--blast_cache false` keeps them in the task instead), so a
  read which leaked in several methods, or in an earlier run, is only searched once. `--store_dir` should be an absolute path
  on a filesystem every task can see;
- a read whose search still fails is written as `failed` and gets `blast_failed = true` in the leak tables. It is counted in
  `num_microbial_map_host_blast_failed_<classifier>` and is neither verified nor shown to be non-human. The run does not stop and
  the summary is still written.
- if a whole BLAST chunk task fails (killed, out of memory, out of time) the reads of that chunk are written as `failed` too,
  and the other chunks of the sample are kept. Such a sample is merged once the other BLAST tasks have finished, which can
  be late in the run.
- the reads that BLAST is asked about (unique, before any `--blast_max_reads` cap) are published as
  `<outdir>/<sample>/<sample>_<method>_blast_queries.fasta.gz`, so the BLAST step can be run again later, for example against a
  local database or after remote searches failed. Read ids match the classifier's, so the hits can be joined back to the leak tables.

### Reads kept for a closer look

With `--publish_reads` (the default) the reads each method calls microbial are copied to
`<outdir>/<sample>/reads/<sample>_<method>_microbial[_1|_2].fastq.gz`. charon does not extract the reads it leaves unclassified,
so for charon these are taken from the input by name (the rows with status `U` in its report) and written as
`<sample>_<method>_unclassified*.fastq.gz`. With `--map_unclassified` (the default) they are also mapped to the reference
(`<sample>/mapped/<sample>.<method>.unclassified.sam`) and summarised per read in `<sample>_<method>_unclassified_alignments.csv.gz`.
These extra steps ignore failures, so they cannot stop a sample.

The mapping columns of the full and leak tables also say where the aligned part of a read is (`query_start`, `query_end`,
in the read as sequenced), `mapq`, `strand`, how many supplementary alignments the read has (`n_supplementary`) and how much of the
read they cover (`supp_aligned_length`), and `total_aligned_prop`, the fraction of the read covered by the primary and
supplementary alignments together. A long read with only a short human piece is then visible as such (chimera or repeat).

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
which are not called host are kept, as for deacon and nohuman, and it is named `charon_host_only` in all outputs. It uses the same charon image
as the normal run, set with `--charon_container`.
The host-only index should be built from the same host sequences and k-mer settings as the full index, so that the only difference
between the two is the microbial category.
