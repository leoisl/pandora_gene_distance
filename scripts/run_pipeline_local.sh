#!/usr/bin/env bash
set -eux

LOG_DIR=logs/
mkdir -p $LOG_DIR

snakemake --use-singularity "$@" >"$LOG_DIR"/pandora_gene_distance.out 2>"$LOG_DIR"/pandora_gene_distance.err

exit 0
