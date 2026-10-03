#!/bin/bash
#$ -l h_rt=12:00:00
#$ -cwd
#$ -q eddie
#$ -pe sharedmem 12
#$ -l h_vmem=4G
#$ -N modkit_medaka

set -euo pipefail

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate modkit202607

BASE=/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/methylation_selected
ASMBASE=/home/s1887840/mazeri_grp/Ana/assemblies/methyl_assemblies
THREADS=${NSLOTS:-12}

for BC in barcode03 barcode14 barcode36 barcode47; do
  ASM="$ASMBASE/medaka_${BC}/consensus_${BC}.fasta"

  for MOD in 6mA 4mC_5mC; do
    BAM="$BASE/mapping_medaka/mapped_bams/$MOD/${BC}.${MOD}.medaka.mapped.sorted.bam"
    OUTDIR="$BASE/modkit_medaka/pileups/$MOD"
    PILEUP="$OUTDIR/${BC}.${MOD}.pileup.bed"
    LOG="$OUTDIR/${BC}.${MOD}.pileup.log"

    mkdir -p "$OUTDIR"

    modkit pileup \
      "$BAM" \
      "$PILEUP" \
      --reference "$ASM" \
      --threads "$THREADS" \
      --log-filepath "$LOG"
  done
done
