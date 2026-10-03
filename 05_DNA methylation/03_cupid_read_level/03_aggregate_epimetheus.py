#!/usr/bin/env python3
import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("input")
p.add_argument("output")
p.add_argument("--read-col", default="read_id")
p.add_argument("--motif-col", default="motif_mod_position")
p.add_argument("--quality-col", default="quality")
args = p.parse_args()

df = pd.read_csv(args.input, sep="\t")
df["methylation_probability"] = pd.to_numeric(
    df[args.quality_col], errors="coerce"
) / 256.0

out = (
    df.groupby(
        [args.read_col, args.motif_col],
        as_index=False,
        dropna=False
    )["methylation_probability"]
    .mean()
)

out.to_csv(args.output, sep="\t", index=False)
