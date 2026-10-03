#!/usr/bin/env python3
import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--amrfinder", required=True)
p.add_argument("--methylation", required=True, nargs="+")
p.add_argument("--output", required=True)
p.add_argument("--read-col-amr", default="Contig id")
p.add_argument("--read-col-methyl", default="read_id")
args = p.parse_args()

amr = pd.read_csv(args.amrfinder, sep="\t")
arg_ids = set(amr[args.read_col_amr].astype(str))

parts = []
for path in args.methylation:
    d = pd.read_csv(path, sep="\t")
    d = d[d[args.read_col_methyl].astype(str).isin(arg_ids)].copy()
    d["methylation_source"] = path
    parts.append(d)

out = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
out.to_csv(args.output, sep="\t", index=False)
