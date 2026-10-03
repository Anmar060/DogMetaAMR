#!/bin/bash
#$ -l h_rt=24:00:00
#$ -cwd
#$ -q eddie
#$ -pe sharedmem 8
#$ -l h_vmem=8G
#$ -N cupid_epi

set -euo pipefail

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate cupid_epimetheus

BC=${1:?barcode}
MOD=${2:?6mA or 4mC_5mC}
MOTIFS=${3:?motif file prepared from NanoMotif bin-motifs.tsv}

BASE=/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/methylation_selected
BAM="$BASE/final_bams/$MOD/${BC}.${MOD}.bam"
OUTDIR="$BASE/cupid_pilot/$BC/read_methylation/$MOD"
FASTQ="$OUTDIR/${BC}.${MOD}.MMML.fastq"
RAW="$OUTDIR/${BC}.${MOD}.epimetheus.tsv"

EPIMETHEUS=/home/s1887840/mazeri_grp/Ana/tools/epimetheus/epimetheus
THREADS=${NSLOTS:-8}

mkdir -p "$OUTDIR"

samtools fastq -T MM,ML "$BAM" > "$FASTQ"

"$EPIMETHEUS" methylation-pattern read-fastq \
  --input "$FASTQ" \
  --output "$RAW" \
  --motifs "$MOTIFS" \
  --threads "$THREADS"
