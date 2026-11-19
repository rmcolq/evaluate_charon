# Analysis scripts

Run these after the pipeline, on its `--outdir`. They need Python with `pandas`, `numpy`, `matplotlib`, `scipy` and
`duckdb` (the `dehumanize` conda env has them). Pool one read type per run (nanopore or illumina).

| Script | What it does |
| --- | --- |
| `leak_regions.py` | Positions of microbial-called reads that map to the human reference (from `intermediate/*.microbial.sam`): per-read, per-locus, per-bin tables and a summary. |
| `compare_common.py` | Host / microbial / unclassified calls on the reads every classifier processed, with 95% Wilson intervals. Reports reads Charon filtered out (below its min quality / length) separately. |
| `plot_leak_ideogram.py` | One figure per chromosome (chr1-22, X, Y, M) and a PDF: repeat-content and gene tracks, CenSat/rDNA and segdups, then a read-density heatmap per classifier. Also writes `enrichment.tsv`. |
| `plot_leak_by_chromosome.py` | Leaked reads per chromosome for each classifier, raw and per 100 Mb, with a chi-square test against length-proportional spread. |

```
scripts/leak_regions.py --outdir out --samples S1 S2 --classifiers charon deacon deacon_sensitive nohuman
scripts/compare_common.py --outdir out --samples S1 S2 --classifiers charon deacon deacon_sensitive nohuman
scripts/plot_leak_ideogram.py --reads results/leak_regions/leak_reads.tsv --samples S1 S2 \
    --classifiers charon deacon deacon_sensitive nohuman --annotations resources/annotations --outdir results/ideogram
scripts/plot_leak_by_chromosome.py --reads results/leak_regions/leak_reads.tsv --samples S1 S2 \
    --classifiers charon deacon deacon_sensitive nohuman -o results/by_chromosome.png
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
