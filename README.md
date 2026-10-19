# evaluate_charon
Nextflow pipeline to run charon and collect statistics


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
