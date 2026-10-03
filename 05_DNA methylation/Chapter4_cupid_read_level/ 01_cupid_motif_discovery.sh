#!/bin/bash
#$ -l h_rt=12:00:00
#$ -cwd
#$ -q eddie
#$ -pe sharedmem 4
#$ -l h_vmem=32G
#$ -N cupid_nanomo

set -euo pipefail

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate modkit202607

BC=${1:?barcode}
MOD=${2:?6mA or 4mC_5mC}

BASE=/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/methylation_selected
ASM=/home/s1887840/mazeri_grp/Ana/assemblies/methyl_assemblies/medaka_${BC}/consensus_${BC}.fasta
PILEUP="$BASE/modkit_medaka/pileups/$MOD/${BC}.${MOD}.pileup.bed"
OUTROOT="$BASE/cupid_pilot/$BC"
MAP="$OUTROOT/inputs/${BC}.contig_as_bin.tsv"
OUT="$OUTROOT/nanomotif_contig/$MOD"

mkdir -p "$(dirname "$MAP")" "$OUT"

python - "$ASM" "$MAP" <<'PY'
from pathlib import Path
import sys
assembly = Path(sys.argv[1])
out = Path(sys.argv[2])
with assembly.open() as fh, out.open("w") as oh:
    for line in fh:
        if line.startswith(">"):
            contig = line[1:].strip().split()[0]
            oh.write(f"{contig}\t{contig}\n")
PY

nanomotif motif_discovery \
  "$ASM" \
  "$PILEUP" \
  -c "$MAP" \
  --out "$OUT" \
  --threads "${NSLOTS:-4}" \
  --threshold_valid_coverage 1 \
  --min_motifs_bin 1 \
  --min_motif_score 0.5 \
  --seed 1 \
  --verbose \
  > "$OUT/run.log" 2>&1
