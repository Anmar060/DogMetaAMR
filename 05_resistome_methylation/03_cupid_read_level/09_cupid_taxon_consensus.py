#!/usr/bin/env python3
import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("top25")
p.add_argument("output")
p.add_argument("--taxid-col", default="host_taxid")
args = p.parse_args()

df = pd.read_csv(args.top25, sep="\t")

out = (
    df.groupby(
        ["query_read_id", args.taxid_col],
        dropna=False
    )
    .agg(
        n_top25_reads=("host_read_id", "size"),
        best_rank=("rank", "min"),
        mean_rank=("rank", "mean"),
        best_rss=("rss", "max"),
        mean_rss=("rss", "mean"),
        mean_rmsd=("rmsd", "mean"),
        mean_shared_motifs=("n_shared_motifs", "mean"),
    )
    .reset_index()
)

out["fraction_top25"] = out["n_top25_reads"] / 25.0
out.to_csv(args.output, sep="\t", index=False)
