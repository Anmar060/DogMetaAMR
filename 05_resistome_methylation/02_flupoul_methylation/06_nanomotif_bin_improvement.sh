#!/bin/bash
#$ -l h_rt=12:00:00
#$ -cwd
#$ -q eddie
#$ -pe sharedmem 8
#$ -l h_vmem=8G
#$ -N nanomotif_improve

set -euo pipefail

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate modkit202607

THREADS=${NSLOTS:-8}
TASK_ID="${SGE_TASK_ID:-}"

if [[ ! "$TASK_ID" =~ ^[1-4]$ ]]; then
  echo "Submit with: qsub -t 1-4 $0" >&2
  exit 1
fi

BARCODES=(barcode03 barcode14 barcode36 barcode47)
BC="${BARCODES[$((TASK_ID - 1))]}"

BASE=/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/methylation_selected
ASMBASE=/home/s1887840/mazeri_grp/Ana/assemblies/methyl_assemblies
ASM="$ASMBASE/medaka_${BC}/consensus_${BC}.fasta"
CONTIG_BIN="$BASE/nanomotif_inputs/contig_bins/${BC}.contig_bin.tsv"

for MOD in 6mA 4mC_5mC; do
  PILEUP="$BASE/modkit_medaka/pileups/$MOD/${BC}.${MOD}.pileup.bed"
  MOTIFS="$BASE/nanomotif_final/$BC/$MOD/flupoul_like/bin-motifs.tsv"
  OUT="$BASE/nanomotif_final/$BC/$MOD/flupoul_like_bin_improvement"

  [[ -s "$MOTIFS" ]] || continue
  mkdir -p "$OUT"

  nanomotif include_contigs \
    --pileup "$PILEUP" \
    --assembly "$ASM" \
    --bin_motifs "$MOTIFS" \
    --contig_bins "$CONTIG_BIN" \
    --run_detect_contamination \
    --threads "$THREADS" \
    --out "$OUT" \
    > "$OUT/run.log" 2>&1
done
