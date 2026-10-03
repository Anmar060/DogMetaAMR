#!/bin/bash
#$ -l h_rt=12:00:00
#$ -cwd
#$ -q eddie
#$ -pe sharedmem 12
#$ -l h_vmem=4G
#$ -N map_modbams_medaka

set -euo pipefail

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate samtools_ana202506

BASE=/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/methylation_selected
ASMBASE=/home/s1887840/mazeri_grp/Ana/assemblies/methyl_assemblies
THREADS=${NSLOTS:-12}

for BC in barcode03 barcode14 barcode36 barcode47; do
  ASM="$ASMBASE/medaka_${BC}/consensus_${BC}.fasta"

  for MOD in 6mA 4mC_5mC; do
    IN="$BASE/final_bams/$MOD/${BC}.${MOD}.bam"
    OUTDIR="$BASE/mapping_medaka/mapped_bams/$MOD"
    OUT="$OUTDIR/${BC}.${MOD}.medaka.mapped.sorted.bam"

    mkdir -p "$OUTDIR"

    samtools fastq -T MM,ML "$IN" \
      | minimap2 -ax map-ont -y -I 16G -t "$THREADS" "$ASM" - \
      | samtools sort -@ "$THREADS" -o "$OUT"

    samtools index -@ "$THREADS" "$OUT"
    samtools flagstat "$OUT" > "${OUT%.bam}.flagstat.txt"
  done
done
