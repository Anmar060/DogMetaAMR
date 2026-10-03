#!/usr/bin/env python3
import argparse
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--kraken", required=True)
p.add_argument("--arg-ids", required=True)
p.add_argument("--output", required=True)
args = p.parse_args()

ids = {x.strip() for x in open(args.arg_ids) if x.strip()}
rows = []

with open(args.kraken) as fh:
    for line in fh:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 3:
            continue
        status, read_id, taxid = fields[:3]
        if read_id in ids:
            rows.append({
                "status": status,
                "read_id": read_id,
                "taxid": taxid,
                "raw": line.rstrip("\n")
            })

pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
