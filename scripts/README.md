# Analysis scripts

Run these after the pipeline, on its `--outdir`, in the conda environment from `environment.yml`. Pool one read type
per run (nanopore or illumina).

| Script | What it does |
| --- | --- |
| `run_hall.sh` | Runs the pipeline on the Hall et al. datasets (`manifests/hall_ont.csv`, `manifests/hall_illumina.csv`), as in our analysis. |
| `leak_regions.py` | Positions of microbial-called reads that map to the human reference (from `<sample>/mapped/*.microbial.sam`): per-read, per-locus, per-bin tables and a summary. |
| `compare_common.py` | Host / microbial / unclassified calls on the reads every classifier processed, with 95% Wilson intervals. Reports reads Charon filtered out (below its min quality / length) separately. |
| `plot_leak_ideogram.py` | One figure per chromosome (chr1-22, X, Y, M) and a PDF: repeat-content and gene tracks, CenSat/rDNA and segdups, then a read-density heatmap per classifier. Also writes `enrichment.tsv`. |
| `charon_param_sweep.py` | For simulated reads with known origin: re-calls a Charon dehost report under other thresholds (a port of `call_host`) and scores each against truth from the Badread read headers. One-at-a-time sweeps, a trade-off plot, a `--confidence` × `--min_proportion_diff` grid and per-feature distributions. Lowering `--min_length`/`--min_quality`, and model parameters, need a rerun. |
| `compare_param_sweeps.py` | Combines the `combination_grid.tsv` of several `charon_param_sweep.py` runs (one per dataset) into one table and a per-dataset plot of human leaks against microbial reads kept. |
| `plot_leak_by_chromosome.py` | Leaked reads per chromosome for each classifier, raw and per 100 Mb, with a chi-square test against length-proportional spread. |

```
scripts/leak_regions.py --outdir out --samples S1 S2 --classifiers charon deacon deacon_sensitive nohuman
scripts/compare_common.py --outdir out --samples S1 S2 --classifiers charon deacon deacon_sensitive nohuman
scripts/plot_leak_ideogram.py --reads results/leak_regions/leak_reads.tsv --samples S1 S2 \
    --classifiers charon deacon deacon_sensitive nohuman --annotations resources/annotations --outdir results/ideogram
scripts/plot_leak_by_chromosome.py --reads results/leak_regions/leak_reads.tsv --samples S1 S2 \
    --classifiers charon deacon deacon_sensitive nohuman -o results/by_chromosome.png
```

`charon_param_sweep.py` runs on the per-read report and log of one `charon dehost` run, not on the pipeline outdir.
Truth is either the simulated FASTQ headers (Badread, or an accession in the read id) or a read_id<TAB>label table:

```
zcat metagenome.ont.fq.gz | awk 'NR%4==1' > metagenome.headers
scripts/charon_param_sweep.py --headers metagenome.headers --report charon_results.tsv --log charon_dehost.log --outdir param_sweep
scripts/charon_param_sweep.py --truth_tsv read2truth.ont.tsv --report charon_results.tsv --log charon_dehost.log --outdir param_sweep
scripts/compare_param_sweeps.py --sweep sim_ont=sim/param_sweep real_ont=real/param_sweep --highlight 7,0.07,0.04,0,0 --outdir comparison
```

`--ref_bed` (default `resources/T2T_ref.bed`) gives the chromosome lengths. The human mitochondrion is the contig
`JAGYVL020000058.1` in the default references and is plotted as chrM.

## Annotation files (not in the repository)

`plot_leak_ideogram.py` reads these from `resources/annotations/`. All come from the UCSC hub for T2T-CHM13v2.0,
`https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/`:

| File | Path under the hub URL |
| --- | --- |
| `GCA_009914755.4.chromAlias.txt` | `GCA_009914755.4.chromAlias.txt` |
| `t2t_censat_CHM13v2.0_trackv2.0.10col.bed` | `bbi/GCA_009914755.4_T2T-CHM13v2.0.censat/` |
| `chm13v1.1.rdna_model.bed` | `bbi/GCA_009914755.4_T2T-CHM13v2.0.rdnaModel/` |
| `T2T-CHM13v2.SDs.bed.gz` | `bbi/GCA_009914755.4_T2T-CHM13v2.0.sedefSegDups/` |
| `chm13v2.0_rmsk.out.gz` (RepeatMasker, 177 MB) | `bbi/GCA_009914755.4_T2T-CHM13v2.0.t2tRepeatMasker/` |
| `GCF_009914755.1_T2T-CHM13v2.0.110.20220412.ncbiRefSeq.gtf.gz` | `genes/` |
