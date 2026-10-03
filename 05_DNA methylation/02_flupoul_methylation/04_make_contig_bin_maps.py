#!/usr/bin/env python3
from pathlib import Path
import re

BIN_ROOT = Path("/home/s1887840/mazeri_grp/Ana/assemblies/methyl_assemblies/metabat2")
OUT_ROOT = Path("/exports/cmvm/eddie/eb/groups/mazeri_grp/Ana/methylation_selected/nanomotif_inputs/contig_bins")
OUT_ROOT.mkdir(parents=True, exist_ok=True)

for barcode in ("barcode03", "barcode14", "barcode36", "barcode47"):
    bindir = BIN_ROOT / f"medaka_metabat2_{barcode}"
    out = OUT_ROOT / f"{barcode}.contig_bin.tsv"
    assignments = {}

    for fasta in sorted(bindir.glob("bin.*.fa")):
        m = re.fullmatch(r"(bin\.\d+)\.fa", fasta.name)
        if not m:
            continue
        bin_id = m.group(1)
        with fasta.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    contig = line[1:].strip().split()[0]
                    if contig in assignments and assignments[contig] != bin_id:
                        raise RuntimeError(f"{contig}: duplicate bin assignment")
                    assignments[contig] = bin_id

    with out.open("w") as handle:
        for contig, bin_id in sorted(assignments.items()):
            handle.write(f"{contig}\t{bin_id}\n")
