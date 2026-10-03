#!/bin/bash
#$ -l h_rt=24:00:00
#$ -cwd
#$ -q eddie
#$ -pe sharedmem 8
#$ -l h_vmem=8G
#$ -N cupid_amrfinder

set -euo pipefail

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate cupid_amrfinder

BC=${1:?barcode}
BAM=${2:?unaligned barcode BAM used for read screening}

BASE=/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/methylation_selected
OUTDIR="$BASE/cupid_pilot/$BC/amrfinder_reads"
FASTQ="$OUTDIR/${BC}.fastq"
FASTA="$OUTDIR/${BC}.fasta"
OUT="$OUTDIR/${BC}.amrfinder.tsv"

mkdir -p "$OUTDIR"

samtools fastq "$BAM" > "$FASTQ"

awk '
  NR % 4 == 1 {sub(/^@/,">"); print}
  NR % 4 == 2 {print}
' "$FASTQ" > "$FASTA"

amrfinder \
  -n "$FASTA" \
  --threads "${NSLOTS:-8}" \
  -o "$OUT"

awk 'NR>1 {print $1}' "$OUT" | sort -u \
  > "$OUTDIR/${BC}.ARG_read_ids.txt"
