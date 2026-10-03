#!/bin/bash
#$ -l h_rt=48:00:00
#$ -cwd
#$ -q eddie
#$ -pe sharedmem 8
#$ -l h_vmem=8G
#$ -N cupid_mob

set -euo pipefail

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate cupid_mobsuite

BC=${1:?barcode}
FASTA=${2:?all-read FASTA}

BASE=/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/methylation_selected
OUTDIR="$BASE/cupid_pilot/$BC/mobsuite_all_reads"

mkdir -p "$OUTDIR"

mob_recon \
  -i "$FASTA" \
  -o "$OUTDIR"
