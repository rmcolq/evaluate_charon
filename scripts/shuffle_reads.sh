#!/usr/bin/env bash
# Shuffle a FASTQ file, or a pair, into a fixed random order. Charon fits its model on the first reads, so a file
# grouped by source genome (as the Hall et al. mixtures are) must be shuffled first. The same seed shuffles both
# files of a pair the same way, so the mates stay in step. The file is held in memory (a few times its gzipped size).
#   scripts/shuffle_reads.sh in.fq.gz out.fq.gz
#   scripts/shuffle_reads.sh in_1.fq.gz in_2.fq.gz out_1.fq.gz out_2.fq.gz
set -euo pipefail
seed=11
if [ $# -eq 2 ]; then
  seqkit shuffle -s "$seed" -o "$2" "$1"
elif [ $# -eq 4 ]; then
  seqkit shuffle -s "$seed" -o "$3" "$1"
  seqkit shuffle -s "$seed" -o "$4" "$2"
else
  echo "usage: $0 in.fq.gz out.fq.gz | in_1.fq.gz in_2.fq.gz out_1.fq.gz out_2.fq.gz" >&2
  exit 1
fi
