# evaluate_charon

A Nextflow pipeline that compares host (human) read removal tools: [Charon](https://github.com/rmcolq/charon), deacon and
nohuman. For each sample it records what each tool calls host and microbial, maps both sets to a human + microbial
reference to find mistakes, and confirms human reads left among the microbial calls ("leaks") with BLAST.

## Setup

1. **Software.** Nextflow 25.10 or later and a container engine (Docker, Singularity or Apptainer). The conda environment
   has Nextflow, Java and the Python packages for the analysis scripts:

   ```bash
   conda env create -f environment.yml
   conda activate evaluate_charon
   ```

   With `-profile conda` the pipeline's processes use this environment too, except Charon, which needs a container.

2. **Reference.** `--refs` is the human (T2T) + EBV + fungi + Mycoplasma reference used for mapping, with `--ref_bed`
   marking its human contigs (`resources/T2T_ref.bed`). The reference will be on Zenodo; until then give its local path.

3. **Charon index.** `--db`: one of the prebuilt indexes on Zenodo, or your own built with `charon index`. The
   prebuilt ones were built with Charon v1.1.7 (from charon-index-build v1.0.0) and need v1.1.7 or later, the
   default `--charon_container`:

   | Read type | Zenodo record | Indexes |
   |---|---|---|
   | Nanopore | [10.5281/zenodo.24559698](https://doi.org/10.5281/zenodo.24559698) | `charon_big.nanopore.idx.gz`, `charon_big_host_only.nanopore.idx.gz`, `charon_small.nanopore.idx.gz` |
   | Illumina | [10.5281/zenodo.24531993](https://doi.org/10.5281/zenodo.24531993) | `charon_big.illumina.idx.gz`, `charon_big_host_only.illumina.idx.gz`, `charon_small.illumina.idx.gz` |

   Give `--db` (and `--db_host_only`) the download URL, `https://zenodo.org/records/<record>/files/<index>?download=1`,
   or a local copy. A gzipped index is decompressed once into `--store_dir`. The deacon and nohuman indexes are
   downloaded by the pipeline on first use.

4. **BLAST database** (for leak confirmation). A local `core_nt` with the taxonomy files, for example:

   ```bash
   mkdir -p /db/blast && cd /db/blast
   update_blastdb.pl --decompress core_nt taxdb
   ```

   Without `--blast_db`, BLAST runs remotely at NCBI, which is slow and sends your reads to NCBI.

Any path option (`--refs`, `--db`, a CSV's reads) can also be a URL: it is downloaded once into `--store_dir`.
For example, `--db https://zenodo.org/records/24559698/files/charon_big.nanopore.idx.gz?download=1`.

## Quick run

```bash
nextflow run main.nf -profile docker \
  --fastq reads.fq.gz --unique_id my_sample --read_type ont \
  --db charon.idx --refs T2T_ref_plus_EBV_fungi_mycoplasma.fa.gz \
  --charon true --deacon true --nohuman true \
  --evaluate_host true --evaluate_microbial true \
  --blast_db /db/blast/core_nt --store_dir $PWD/store_dir --outdir results
```

Input, one of:

| Option | Input |
| --- | --- |
| `--fastq` (and `--fastq2` for a pair) | one sample; `--unique_id` names it |
| `--fastq_csv` | `sample_id,fastq_1,fastq_2` (`fastq_2` empty for single-end reads), plus an optional `shuffle` column: `true` shuffles that sample's reads first |
| `--bam_csv` | `sample_id,filepath` of BAM/CRAM files; add `--paired true` for read pairs and `--cram_ref` for CRAMs |
| `--fastq_dir` | one single-end FASTQ per sample |

Main options:

| Option | Default | |
| --- | --- | --- |
| `--charon`, `--deacon`, `--deacon_sensitive`, `--nohuman`, `--charon_host_only` | `false` | the tools to compare |
| `--read_type` | `ont` | `ont` or `illumina` (pairs are always mapped as short reads) |
| `--charon_args` | | extra options for `charon dehost`, e.g. `"--confidence 10"` |
| `--evaluate_host`, `--evaluate_microbial` | `false` | map the host and microbial calls to find mistakes |
| `--blast_db` | remote | BLAST database path without extensions |
| `--max_concurrent_samples` | `1` | samples in flight at once, to limit disk use |
| `--max_cpus`, `--max_memory` | `16`, `128.GB` | caps for any one task; `-profile local` suits a laptop |

## Outputs

Per sample, in `<outdir>/<sample>/`:

- `<sample>_<tool>_summary.csv`: counts and rates for each tool: host, microbial and unclassified calls, host reads that
  map elsewhere, and leaks by tier;
- `<sample>_<tool>_leaks.csv.gz`: every read called microbial, with its features, mapping, BLAST hits and leak tier
  (1 called microbial, 2 maps to human, 3 BLAST human hit, 4 human-specific BLAST hit);
- `<sample>_<tool>_calls.csv.gz`: every read and its call (both are `.parquet` instead when `pyarrow` is installed);
- `blast/`, `mapped/`, `reads/`: BLAST results, alignments and the reads each tool kept.

`<outdir>/pipeline_info/` has the run report and the time and memory of every task. See [docs/details.md](docs/details.md)
for the leak tiers, BLAST, paired reads and other details.

## Replicating our analyses

**Hall et al. metagenomes.** The artificial and simulated metagenomes from Hall et al., nanopore and Illumina, are in
`manifests/hall_ont.csv` and `manifests/hall_illumina.csv` (artificial and simulated in each), which point at the original reads on Zenodo. The originals are grouped by source genome and
Charon fits its model on the first reads, so these manifests set `shuffle=true`: the pipeline shuffles each file
(`seqkit shuffle`, seqkit 2.14.0, seed `--shuffle_seed`, default 11; the same seed for both files of a pair keeps the
mates in step) before running the classifiers. Shuffling holds the file in memory, up to a few times its gzipped size.
`scripts/run_hall.sh` runs both:

```bash
BLAST_DB=/db/blast/core_nt scripts/run_hall.sh results/hall
```

By default it uses the `charon_big` and `charon_big_host_only` indexes from Zenodo; set `CHARON_DB_ONT`,
`CHARON_DB_ILLUMINA`, `CHARON_DB_HOST_ONLY_ONT` and `CHARON_DB_HOST_ONLY_ILLUMINA` to use others (local paths or URLs).

The reads come from these Zenodo records:

| Dataset | Zenodo record | Files |
|---|---|---|
| Simulated, Illumina | [8339791](https://zenodo.org/records/8339791) | `metagenome_R1.illumina.fq.gz`, `metagenome_R2.illumina.fq.gz` |
| Simulated, nanopore | [8339789](https://zenodo.org/records/8339789) | `metagenome.ont.fq.gz` |
| Artificial, both | [10472796](https://zenodo.org/records/10472796) | `all_1.illumina.fq.gz`, `all_2.illumina.fq.gz`, `all.ont.fq.gz`, `read2truth.{illumina,ont}.tsv` |

Add `shuffle=true` to your own `--fastq_csv` rows to shuffle other
read files the same way.

**Human samples (HPRC).** `manifests/human_nanopore.csv` (and a 10-sample `human_nanopore_mini.csv`) and
`manifests/human_illumina_paired.csv` list public HPRC reads on S3. They are HPRC release 2 samples chosen so that none of
their assemblies is in the Charon index: the indexes are built from the
[charon-index-build](https://github.com/rmcolq/charon-index-build) manifests (`manifests/<profile>/human_assemblies.tsv`),
which leave these samples out. Reads from a sample whose assembly is in the index would test recall of sequence the
index already holds, not of unseen human genomes.

```bash
nextflow run main.nf -profile docker --bam_csv manifests/human_nanopore_mini.csv --s3_anonymous true \
  --read_type ont --db "https://zenodo.org/records/24559698/files/charon_big.nanopore.idx.gz?download=1" --charon true --deacon true --deacon_sensitive true --nohuman true \
  --evaluate_host true --evaluate_microbial true --blast_db /db/blast/core_nt --store_dir $PWD/store_dir --outdir results/hprc_ont

nextflow run main.nf -profile docker --bam_csv manifests/human_illumina_paired.csv --paired true --s3_anonymous true \
  --cram_ref GRCh38_full_analysis_set_plus_decoy_hla.fa --read_type illumina \
  --db "https://zenodo.org/records/24531993/files/charon_big.illumina.idx.gz?download=1" \
  --charon true --deacon true --deacon_sensitive true --nohuman true \
  --evaluate_host true --evaluate_microbial true --blast_db /db/blast/core_nt --store_dir $PWD/store_dir --outdir results/hprc_illumina
```

The Illumina CRAMs are CRAM 3.0 without an embedded reference, so they need the 1000 Genomes GRCh38 reference they
were aligned to (`GRCh38_full_analysis_set_plus_decoy_hla.fa`, named in their `@SQ UR` tags).

**Figures and tables.** The analysis scripts in `scripts/` run on the pipeline's output; see
[scripts/README.md](scripts/README.md).
