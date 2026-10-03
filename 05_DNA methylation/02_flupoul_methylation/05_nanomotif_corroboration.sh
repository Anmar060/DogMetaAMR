#!/bin/bash
#$ -l h_rt=12:00:00
#$ -cwd
#$ -q eddie
#$ -pe sharedmem 8
#$ -l h_vmem=8G
#$ -N nanomotif_compare

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
OUTBASE="$BASE/nanomotif_final/$BC"

run_nm () {
  local setting=$1
  local mod=$2
  local pileup=$3
  local cov=$4
  local score=$5
  local out="$OUTBASE/$mod/$setting"

  mkdir -p "$out"

  nanomotif motif_discovery \
    "$ASM" \
    "$pileup" \
    -c "$CONTIG_BIN" \
    --out "$out" \
    --threads "$THREADS" \
    --threshold_valid_coverage "$cov" \
    --min_motif_score "$score" \
    --min_motifs_bin 50 \
    --methylation_threshold_low 0.3 \
    --methylation_threshold_high 0.7 \
    --seed 1 \
    --verbose \
    > "$out/run.log" 2>&1
}

for MOD in 6mA 4mC_5mC; do
  run_nm \
    flupoul_like "$MOD" \
    "$BASE/modkit_medaka/pileups/$MOD/${BC}.${MOD}.pileup.bed" \
    1 0.2

  run_nm \
    conservative "$MOD" \
    "$BASE/modkit_medaka_nanofilter/pileups/$MOD/${BC}.${MOD}.nanofilter.pileup.bed" \
    5 1.5
done
