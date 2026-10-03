#!/usr/bin/env python3
"""
rank_gene_variability.py - rank genes by how much their coverage varies across samples,
and (optionally) test whether that variation lines up with sequencing batch.

Input is the long table written by summarize_depth.sh when bam_depth.sh was run with --gene-bed:
    sample_id  gene  ensg  n_bases  mean_depth  [frac_10x  frac_20x ...]

Usage:
    rank_gene_variability.py gene_depth.tsv [--batch-map BATCHES.tsv] [--outdir DIR]
                             [--metric mean_depth] [--rank-by norm_cv] [--no-normalize]

Outputs (in --outdir, default: next to the input):
    gene_variability.tsv     every gene, ranked most -> least variable across samples
    gene_batch_effects.tsv   (with --batch-map) per-gene test of coverage vs batch, ranked
    gene_batch_means.tsv     (with --batch-map) mean coverage of every gene in every batch

Batch map: two columns, sample_id and batch (tab/comma/space separated; a header line is optional).
sample_id must match the sample_id column of the depth table (BAM file name without .bam).

Requires python3 with numpy, pandas, scipy.
"""

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)  # all-NaN rows etc. are handled explicitly

RANK_CHOICES = ["norm_cv", "cv_excess", "norm_sd", "cv", "sd"]


def log(msg):
    print(msg, file=sys.stderr, flush=True)


def die(msg):
    sys.exit(f"Error: {msg}")


# --------------------------------------------------------------------------- reading

def read_depth_table(path, metric):
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    for col in ("sample_id", "gene", metric):
        if col not in header:
            die(f"column '{col}' not found in {path} (columns: {', '.join(header)})")
    frac_cols = [c for c in header if c.startswith("frac_")]
    cols = ["sample_id", "gene"] + [c for c in ("ensg", "n_bases") if c in header]
    cols += [c for c in dict.fromkeys([metric] + frac_cols) if c not in cols]
    dtypes = {c: "category" for c in ("sample_id", "gene", "ensg") if c in cols}
    dtypes.update({c: "float32" for c in cols if c not in dtypes})
    df = pd.read_csv(path, sep="\t", usecols=cols, dtype=dtypes, keep_default_na=False)
    return df, frac_cols


def build_matrices(df, metric, extra_cols):
    """Return gene info table, sample names, and {column: genes x samples float matrix}."""
    s_codes = df["sample_id"].cat.codes.to_numpy().astype(np.int64)
    samples = np.asarray(df["sample_id"].cat.categories.astype(str))
    g_codes = df["gene"].cat.codes.to_numpy().astype(np.int64)
    if "ensg" in df:
        e_codes = df["ensg"].cat.codes.to_numpy().astype(np.int64)
        combined = g_codes * (int(e_codes.max()) + 1) + e_codes
    else:
        combined = g_codes
    _, first_idx, gi = np.unique(combined, return_index=True, return_inverse=True)
    n_genes, n_samples = len(first_idx), len(samples)

    info = pd.DataFrame({"gene": df["gene"].iloc[first_idx].astype(str).to_numpy()})
    if "ensg" in df:
        info["ensg"] = df["ensg"].iloc[first_idx].astype(str).to_numpy()
    if "n_bases" in df:
        info["n_bases"] = df["n_bases"].to_numpy()[first_idx].astype(np.int64)

    mats = {}
    for col in [metric] + extra_cols:
        M = np.full((n_genes, n_samples), np.nan, dtype=np.float32)
        M[gi, s_codes] = df[col].to_numpy()
        mats[col] = M
    return info, samples, mats


# --------------------------------------------------------------------------- variability

def normalize_by_sample(M):
    """Rescale each sample so its median gene depth equals the overall median of those medians.
    Removes whole-sample depth differences (library size / capture efficiency) so that only
    gene-specific differences remain."""
    pos = np.where(M > 0, M, np.nan)
    size = np.nanmedian(pos, axis=0)
    size = np.where(np.isfinite(size) & (size > 0), size, np.nan)
    return M * (np.nanmedian(size) / size)[None, :], size


def row_summary(M, min_samples):
    n = np.sum(~np.isnan(M), axis=1)
    mean = np.nanmean(M, axis=1)
    sd = np.nanstd(M, axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = np.where(mean > 0, sd / mean, np.nan)
    ok = n >= min_samples
    return n, np.where(ok, mean, np.nan), np.where(ok, sd, np.nan), np.where(ok, cv, np.nan)


def cv_excess(norm_mean, norm_cv, n_bins=20):
    """log(CV) minus the median log(CV) of genes with similar mean depth.
    Low-coverage genes have a high CV from counting noise alone; this asks which genes vary
    more than other genes at the same coverage level."""
    out = np.full(len(norm_cv), np.nan)
    ok = np.isfinite(norm_cv) & (norm_cv > 0) & np.isfinite(norm_mean) & (norm_mean > 0)
    if ok.sum() < 2 * n_bins:
        return out
    lcv = np.log(norm_cv[ok])
    bins = pd.qcut(np.log(norm_mean[ok]), q=n_bins, duplicates="drop")
    med = pd.Series(lcv).groupby(bins, observed=True).transform("median").to_numpy()
    out[ok] = lcv - med
    return out


# --------------------------------------------------------------------------- batch effects

def read_batch_map(path):
    raw = pd.read_csv(path, sep=None, engine="python", header=None, dtype=str, comment="#",
                      skipinitialspace=True)
    if raw.shape[1] < 2:
        die(f"{path}: need two columns (sample_id, batch)")
    raw = raw.iloc[:, :2].apply(lambda c: c.str.strip())
    raw.columns = ["sample_id", "batch"]
    if raw.iloc[0, 0].lower() in ("sample_id", "sample", "id", "sample_name"):
        raw = raw.iloc[1:]
    raw["sample_id"] = raw["sample_id"].str.replace(r"\.bam$", "", regex=True)
    raw = raw.dropna()
    dup = raw["sample_id"][raw["sample_id"].duplicated()]
    if len(dup):
        die(f"{path}: sample(s) listed more than once: {', '.join(dup.unique()[:5])}")
    return dict(zip(raw["sample_id"], raw["batch"]))


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(len(pv))
    res[order] = np.minimum(ranked, 1.0)
    out[ok] = res
    return out


def batch_analysis(info, M, samples, batch_of, min_batch_size, min_samples):
    sample_batch = np.array([batch_of.get(s) for s in samples], dtype=object)
    n_unmapped = int(np.sum([b is None for b in sample_batch]))
    if n_unmapped:
        missing = [s for s, b in zip(samples, sample_batch) if b is None]
        log(f"Warning: {n_unmapped} sample(s) in the depth table are not in the batch map and are "
            f"excluded from the batch analysis (e.g. {', '.join(missing[:5])})")
    in_map = set(samples)
    extra = [s for s in batch_of if s not in in_map]
    if extra:
        log(f"Warning: {len(extra)} sample(s) in the batch map are not in the depth table "
            f"(e.g. {', '.join(extra[:5])})")

    counts = pd.Series([b for b in sample_batch if b is not None]).value_counts()
    keep = [b for b, c in counts.items() if c >= min_batch_size]
    dropped = {b: int(c) for b, c in counts.items() if c < min_batch_size}
    if dropped:
        log(f"Warning: dropping batches with fewer than {min_batch_size} samples: {dropped}")
    if len(keep) < 2:
        die("need at least 2 batches with enough samples for a batch analysis")
    keep = sorted(keep, key=lambda b: (len(str(b)), str(b)))
    log("Batches used: " + ", ".join(f"{b} (n={counts[b]})" for b in keep))

    masks = [(sample_batch == b) for b in keep]
    use = np.any(masks, axis=0)
    X = M[:, use]
    masks = [m[use] for m in masks]
    G = X.shape[0]

    # per-batch means / counts
    n_b = np.stack([np.sum(~np.isnan(X[:, m]), axis=1) for m in masks], axis=1)
    mean_b = np.stack([np.nanmean(X[:, m], axis=1) for m in masks], axis=1)
    n_tot = n_b.sum(axis=1)
    grand = np.nanmean(X, axis=1)
    ss_between = np.nansum(n_b * (mean_b - grand[:, None]) ** 2, axis=1)
    ss_total = np.nansum((X - grand[:, None]) ** 2, axis=1)
    k = len(keep)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta2 = np.where(ss_total > 0, ss_between / ss_total, np.nan)
        eta2_adj = 1 - (1 - eta2) * (n_tot - 1) / (n_tot - k)   # adjusted R^2: removes the chance
        # fit that grows with the number of batches
        fold = np.nanmax(mean_b, axis=1) / np.nanmin(np.where(mean_b > 0, mean_b, np.nan), axis=1)

    # Kruskal-Wallis per gene (rank-based, so robust to outlier samples)
    H = np.full(G, np.nan)
    P = np.full(G, np.nan)
    for i in range(G):
        groups = [X[i, m][~np.isnan(X[i, m])] for m in masks]
        if sum(len(g) >= 1 for g in groups) < 2 or sum(len(g) for g in groups) < min_samples:
            continue
        allv = np.concatenate(groups)
        if np.all(allv == allv[0]):
            H[i], P[i] = 0.0, 1.0
            continue
        try:
            H[i], P[i] = stats.kruskal(*[g for g in groups if len(g)])
        except ValueError:
            pass
    fdr = bh_fdr(P)

    top = np.argmax(np.where(np.isnan(mean_b), -np.inf, mean_b), axis=1)
    bottom = np.argmin(np.where(np.isnan(mean_b) | (mean_b <= 0), np.inf, mean_b), axis=1)
    res = info.copy()
    res["n_samples"] = n_tot
    res["n_batches"] = k
    res["kw_H"] = H
    res["kw_p"] = P
    res["kw_fdr"] = fdr
    res["eta2"] = eta2
    res["eta2_adj"] = eta2_adj
    res["fold_range"] = fold
    res["highest_batch"] = np.array(keep, dtype=object)[top]
    res["lowest_batch"] = np.array(keep, dtype=object)[bottom]
    means = info.copy()
    for j, b in enumerate(keep):
        means[f"mean_{b}"] = mean_b[:, j]
    return res, means, keep, sample_batch


# --------------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gene_depth", help="gene_depth.tsv from summarize_depth.sh")
    ap.add_argument("--batch-map", help="file mapping sample_id to batch (enables the batch analysis)")
    ap.add_argument("--outdir", help="output directory (default: directory of the input file)")
    ap.add_argument("--metric", default="mean_depth",
                    help="column to analyse (default mean_depth; e.g. frac_20x)")
    ap.add_argument("--rank-by", default="norm_cv", choices=RANK_CHOICES,
                    help="statistic used to rank genes (default norm_cv = sd/mean after normalization)")
    ap.add_argument("--no-normalize", action="store_true",
                    help="skip rescaling samples to a common median depth (mean_depth only)")
    ap.add_argument("--min-samples", type=int, default=3,
                    help="genes with fewer samples than this get no statistics (default 3)")
    ap.add_argument("--min-batch-size", type=int, default=3,
                    help="batches with fewer samples than this are dropped (default 3)")
    ap.add_argument("--min-eta2", type=float, default=0.1,
                    help="adjusted eta^2 threshold for flagging a batch effect (default 0.1)")
    ap.add_argument("--fdr", type=float, default=0.05, help="FDR threshold for flagging (default 0.05)")
    args = ap.parse_args()

    if not os.path.isfile(args.gene_depth):
        die(f"not found: {args.gene_depth}")
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.gene_depth))
    os.makedirs(outdir, exist_ok=True)

    log(f"Reading {args.gene_depth} ...")
    df, frac_cols = read_depth_table(args.gene_depth, args.metric)
    extra = [c for c in frac_cols if c != args.metric]
    info, samples, mats = build_matrices(df, args.metric, extra)
    del df
    M = mats[args.metric]
    log(f"{M.shape[0]} genes x {M.shape[1]} samples")
    if M.shape[1] < args.min_samples:
        die("fewer samples than --min-samples")

    # normalization only makes sense for depth, not for fractions
    normalize = (args.metric == "mean_depth") and not args.no_normalize
    if normalize:
        N, size = normalize_by_sample(M)
        log(f"Normalized: each sample rescaled to the median sample's depth "
            f"(median gene depth per sample: min {np.nanmin(size):.1f}, median {np.nanmedian(size):.1f}, "
            f"max {np.nanmax(size):.1f})")
    else:
        N = M
        log("Normalization: not applied" + ("" if args.metric == "mean_depth" else f" (metric is {args.metric})"))

    # ---- variability table
    n, mean, sd, cv = row_summary(M, args.min_samples)
    _, nmean, nsd, ncv = row_summary(N, args.min_samples)
    var = info.copy()
    var["n_samples"] = n
    var["mean"], var["sd"], var["cv"] = mean, sd, cv
    var["norm_mean"], var["norm_sd"], var["norm_cv"] = nmean, nsd, ncv
    var["cv_excess"] = cv_excess(nmean, ncv)
    for c in extra:
        _, fm, fs, _ = row_summary(mats[c], args.min_samples)
        var[f"mean_{c}"], var[f"sd_{c}"] = fm, fs
    key = var[args.rank_by]
    var.insert(0, "rank", key.rank(ascending=False, method="first", na_option="bottom").astype(int))
    var = var.sort_values("rank")
    path = os.path.join(outdir, "gene_variability.tsv")
    var.to_csv(path, sep="\t", index=False, float_format="%.5g", na_rep="NA")
    log(f"Wrote {path}  (ranked by {args.rank_by}, highest first)")
    cols = [c for c in ["rank", "gene", "ensg", "norm_mean", "norm_sd", "norm_cv", "cv_excess"] if c in var]
    log("\nMost variable genes:\n" + var[cols].head(10).to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    # ---- batch analysis
    if args.batch_map:
        if not os.path.isfile(args.batch_map):
            die(f"batch map not found: {args.batch_map}")
        batch_of = read_batch_map(args.batch_map)
        res, means, keep, sample_batch = batch_analysis(info, N, samples, batch_of,
                                                        args.min_batch_size, args.min_samples)
        # whole-sample depth by batch, so a global batch difference is visible even though it is normalized away
        if normalize:
            pb = pd.Series(size, index=samples).groupby(pd.Series(sample_batch, index=samples)).median()
            log("\nMedian per-sample gene depth by batch (removed by normalization):\n"
                + pb.reindex(keep).to_string(float_format=lambda x: f"{x:.1f}"))
        res["flag"] = (res["kw_fdr"] < args.fdr) & (res["eta2_adj"] >= args.min_eta2)
        res = res.sort_values(["eta2_adj"], ascending=False, na_position="last")
        res.insert(0, "rank", np.arange(1, len(res) + 1))
        path = os.path.join(outdir, "gene_batch_effects.tsv")
        res.to_csv(path, sep="\t", index=False, float_format="%.5g", na_rep="NA")
        mpath = os.path.join(outdir, "gene_batch_means.tsv")
        means.loc[res.index].to_csv(mpath, sep="\t", index=False, float_format="%.5g", na_rep="NA")
        log(f"\nWrote {path}\nWrote {mpath}")
        n_sig = int((res["kw_fdr"] < args.fdr).sum())
        n_flag = int(res["flag"].sum())
        log(f"Genes with Kruskal-Wallis FDR < {args.fdr}: {n_sig} of {len(res)}; "
            f"also adjusted eta^2 >= {args.min_eta2}: {n_flag}")
        cols = [c for c in ["rank", "gene", "ensg", "eta2_adj", "kw_fdr", "fold_range",
                            "highest_batch", "lowest_batch"] if c in res]
        log("\nStrongest batch effects:\n" + res[cols].head(10).to_string(index=False, float_format=lambda x: f"{x:.3g}"))


if __name__ == "__main__":
    main()
