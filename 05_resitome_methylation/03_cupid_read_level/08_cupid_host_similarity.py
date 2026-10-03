#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--profiles", required=True)
p.add_argument("--queries", required=True)
p.add_argument("--candidates", required=True)
p.add_argument("--output-all", required=True)
p.add_argument("--output-top25", required=True)
p.add_argument("--read-col", default="read_id")
p.add_argument("--motif-col", default="motif_mod_position")
p.add_argument("--prob-col", default="methylation_probability")
args = p.parse_args()

profiles = pd.read_csv(args.profiles, sep="\t")
queries = pd.read_csv(args.queries, sep="\t")
candidates = pd.read_csv(args.candidates, sep="\t")

profile_map = {
    str(rid): g.set_index(args.motif_col)[args.prob_col].to_dict()
    for rid, g in profiles.groupby(args.read_col)
}

rows = []

for _, q in queries.iterrows():
    qid = str(q[args.read_col])
    qp = profile_map.get(qid, {})
    if not qp:
        continue

    for _, c in candidates.iterrows():
        cid = str(c[args.read_col])
        cp = profile_map.get(cid, {})
        shared = sorted(set(qp) & set(cp))
        if not shared:
            continue

        diff = np.array([qp[m] - cp[m] for m in shared], dtype=float)
        rmsd = float(np.sqrt(np.mean(diff ** 2)))
        rss = max(0.0, 1.0 - rmsd) * len(shared)

        row = {
            "query_read_id": qid,
            "host_read_id": cid,
            "n_shared_motifs": len(shared),
            "rmsd": rmsd,
            "rss": rss
        }

        for col in c.index:
            if col not in row:
                row[col] = c[col]

        rows.append(row)

scores = pd.DataFrame(rows)
scores.to_csv(args.output_all, sep="\t", index=False)

top25 = (
    scores.sort_values(
        ["query_read_id", "rss", "rmsd"],
        ascending=[True, False, True]
    )
    .groupby("query_read_id", group_keys=False)
    .head(25)
    .copy()
)

top25["rank"] = top25.groupby("query_read_id").cumcount() + 1
top25.to_csv(args.output_top25, sep="\t", index=False)
