#!/usr/bin/env bash
# Replicate the Hall et al. runs: the artificial and simulated metagenomes, one run for nanopore and one for Illumina.
#   scripts/run_hall.sh <output directory> [extra nextflow options]
# Set the paths below first. The Charon indexes and the reference can be local paths or URLs (downloaded once into
# --store_dir). BLAST needs a local database: core_nt matches the web BLAST used to curate the truth sets.
set -euo pipefail

CHARON_DB_ONT=${CHARON_DB_ONT:?set to the nanopore charon index, e.g. charon_big.nanopore.idx}
CHARON_DB_ILLUMINA=${CHARON_DB_ILLUMINA:?set to the Illumina charon index, e.g. charon_big.illumina.idx}
BLAST_DB=${BLAST_DB:?set to the BLAST database path without extensions, e.g. /db/blast/core_nt}
STORE_DIR=${STORE_DIR:-$PWD/store_dir}
PROFILE=${PROFILE:-docker}

outdir=${1:?usage: $0 <output directory> [nextflow options]}
shift
pipeline=$(cd "$(dirname "$0")/.." && pwd)

run() {
    local name=$1 manifest=$2 db=$3 read_type=$4
    shift 4
    mkdir -p "$outdir/$name"
    nextflow run "$pipeline/main.nf" -profile "$PROFILE" \
        --fastq_csv "$pipeline/manifests/$manifest" \
        --db "$db" --read_type "$read_type" \
        --charon true --deacon true --deacon_sensitive true --nohuman true \
        --evaluate_host true --evaluate_microbial true \
        --blast_db "$BLAST_DB" --store_dir "$STORE_DIR" \
        --outdir "$outdir/$name" \
        -ansi-log false "$@" > "$outdir/$name/nextflow.log"
}

run ont hall_ont.csv "$CHARON_DB_ONT" ont "$@"
run illumina hall_illumina.csv "$CHARON_DB_ILLUMINA" illumina "$@"
