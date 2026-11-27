#!/usr/bin/env bash
# Replicate the Hall et al. runs: the artificial and simulated metagenomes, one run for nanopore and one for Illumina.
#   scripts/run_hall.sh <output directory> [extra nextflow options]
# The Charon indexes default to the charon_big indexes on Zenodo; they and the reference can be local paths or URLs
# (downloaded once into --store_dir). BLAST needs a local database: core_nt matches the web BLAST used to curate the truth sets.
# Each run has its own --slot_dir, so it does not clash with another pipeline run started from the same directory.
set -euo pipefail

ONT_INDEXES=https://zenodo.org/records/24559698/files        # 10.5281/zenodo.24559698
ILLUMINA_INDEXES=https://zenodo.org/records/24531993/files   # 10.5281/zenodo.24531993
CHARON_DB_ONT=${CHARON_DB_ONT:-$ONT_INDEXES/charon_big.nanopore.idx.gz?download=1}
CHARON_DB_ILLUMINA=${CHARON_DB_ILLUMINA:-$ILLUMINA_INDEXES/charon_big.illumina.idx.gz?download=1}
CHARON_DB_HOST_ONLY_ONT=${CHARON_DB_HOST_ONLY_ONT:-$ONT_INDEXES/charon_big_host_only.nanopore.idx.gz?download=1}
CHARON_DB_HOST_ONLY_ILLUMINA=${CHARON_DB_HOST_ONLY_ILLUMINA:-$ILLUMINA_INDEXES/charon_big_host_only.illumina.idx.gz?download=1}
BLAST_DB=${BLAST_DB:?set to the BLAST database path without extensions, e.g. /db/blast/core_nt}
STORE_DIR=${STORE_DIR:-$PWD/store_dir}
PROFILE=${PROFILE:-docker}

outdir=${1:?usage: $0 <output directory> [nextflow options]}
shift
mkdir -p "$outdir" && outdir=$(cd "$outdir" && pwd)  # absolute, as the slot dir is used from task directories
pipeline=$(cd "$(dirname "$0")/.." && pwd)

run() {
    local name=$1 manifest=$2 db=$3 db_host_only=$4 read_type=$5
    shift 5
    mkdir -p "$outdir/$name"
    nextflow run "$pipeline/main.nf" -profile "$PROFILE" \
        --fastq_csv "$pipeline/manifests/$manifest" \
        --db "$db" --db_host_only "$db_host_only" --read_type "$read_type" \
        --charon true --charon_host_only true --deacon true --deacon_sensitive true --nohuman true \
        --evaluate_host true --evaluate_microbial true \
        --blast_db "$BLAST_DB" --store_dir "$STORE_DIR" \
        --outdir "$outdir/$name" --slot_dir "$outdir/$name/.charon_slots" \
        -ansi-log false "$@" > "$outdir/$name/nextflow.log"
}

run ont hall_ont.csv "$CHARON_DB_ONT" "$CHARON_DB_HOST_ONLY_ONT" ont "$@"
run illumina hall_illumina.csv "$CHARON_DB_ILLUMINA" "$CHARON_DB_HOST_ONLY_ILLUMINA" illumina "$@"
