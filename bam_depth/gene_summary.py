#!/usr/bin/env python3
"""
gene_summary.py - one row per gene, summarized across samples.

Usage: gene_summary.py GENE_DEPTH.tsv OUT.tsv
    GENE_DEPTH.tsv  long table from summarize_depth.sh: sample_id, gene, ensg, n_bases, mean_depth[, frac_*]
    OUT.tsv         columns: gene, ensg, n_bases, n_samples,
                             mean_depth, sd_depth, cv_depth, median_depth, min_depth, max_depth,
                             n_lt_1x, n_lt_10x   (samples whose gene depth is below 1x / 10x),
                             and mean_/sd_ of each frac_* column (fraction of the gene's bases at >= Nx)
    sd is the sample standard deviation (n-1); cv = sd / mean.

Called by summarize_depth.sh; can also be run on its own. Requires python3 with pandas.
"""

import sys

import numpy as np
import pandas as pd


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1:]

    header = pd.read_csv(src, sep="\t", nrows=0).columns.tolist()
    for col in ("gene", "mean_depth"):
        if col not in header:
            sys.exit(f"Error: column '{col}' not found in {src}")
    keys = [c for c in ("gene", "ensg") if c in header]
    fracs = [c for c in header if c.startswith("frac_")]
    cols = keys + [c for c in ("n_bases",) if c in header] + ["mean_depth"] + fracs
    dtypes = {c: "category" for c in keys}
    dtypes.update({c: "float32" for c in cols if c not in keys})
    df = pd.read_csv(src, sep="\t", usecols=cols, dtype=dtypes, keep_default_na=False)
    if df.empty:
        sys.exit(f"Error: no rows in {src}")

    df["lt1"] = df["mean_depth"] < 1
    df["lt10"] = df["mean_depth"] < 10
    g = df.groupby(keys, observed=True, sort=True)

    out = pd.DataFrame(index=g.size().index)
    if "n_bases" in df:
        out["n_bases"] = g["n_bases"].first().astype("int64")
    out["n_samples"] = g["mean_depth"].count()
    out["mean_depth"] = g["mean_depth"].mean()
    out["sd_depth"] = g["mean_depth"].std(ddof=1)
    out["cv_depth"] = out["sd_depth"] / out["mean_depth"].where(out["mean_depth"] > 0)
    out["median_depth"] = g["mean_depth"].median()
    out["min_depth"] = g["mean_depth"].min()
    out["max_depth"] = g["mean_depth"].max()
    out["n_lt_1x"] = g["lt1"].sum().astype("int64")
    out["n_lt_10x"] = g["lt10"].sum().astype("int64")
    for c in fracs:
        out[f"mean_{c}"] = g[c].mean()
        out[f"sd_{c}"] = g[c].std(ddof=1)

    out = out.reset_index()
    out.to_csv(dst, sep="\t", index=False, float_format="%.4f", na_rep="NA")
    print(f"Gene summary written: {dst} ({len(out)} genes)")


if __name__ == "__main__":
    main()
