#!/bin/bash
set -euo pipefail

ASMBASE=/home/s1887840/mazeri_grp/Ana/assemblies/methyl_assemblies
OUTDIR=${1:-.}

mkdir -p "$OUTDIR"

abricate --version > "$OUTDIR/abricate_version.txt"
abricate --list > "$OUTDIR/abricate_database_versions.tsv"

abricate \
  --db megares \
  "$ASMBASE/medaka_barcode03/consensus_barcode03.fasta" \
  "$ASMBASE/medaka_barcode14/consensus_barcode14.fasta" \
  "$ASMBASE/medaka_barcode36/consensus_barcode36.fasta" \
  "$ASMBASE/medaka_barcode47/consensus_barcode47.fasta" \
  > "$OUTDIR/abricate_megares_master.tsv"

abricate --summary "$OUTDIR/abricate_megares_master.tsv" \
  > "$OUTDIR/abricate_megares_summary.tsv"
