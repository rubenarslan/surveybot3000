"""Rank-based recall of related scales, documented-sign pooling (variant B),
cross-instrument pairs, from the existing scale_pair_predictions.csv.

    python pooling_comparison/rank_recall.py

For each query scale, all cross-instrument database scales are ranked by
absolute predicted cosine. For every related pair (|empirical r| >= .40, and
separately >= .30) the rank of the database scale in the query's list and the
list size are recorded (each pair enters twice, once per direction). Reported:
share of related scales within the top 5/10/25/50/100, median rank, and the
chance reference under random ordering (min(k, list size) / list size, averaged
over pairs). Writes rank_recall.csv, rank_recall.md and rank_recall_<dataset>.png
to results/pooling_comparison/. Existing outputs are not modified.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results", "pooling_comparison")
VARIANT = "documented"
THRESHOLDS = [0.40, 0.30]
KS = [5, 10, 25, 50, 100]


def ranks_for(d: pd.DataFrame) -> pd.DataFrame:
    cols = ["query", "db", "empirical_r", VARIANT]
    long = pd.concat([
        d.rename(columns={"scale_a": "query", "scale_b": "db"})[cols],
        d.rename(columns={"scale_b": "query", "scale_a": "db"})[cols],
    ], ignore_index=True)
    long["abs_pred"] = long[VARIANT].abs()
    long["abs_emp"] = long.empirical_r.abs()
    # rank 1 = highest |predicted| within the query's list; ties broken by order
    long["rank"] = long.groupby("query")["abs_pred"].rank(method="first", ascending=False).astype(int)
    long["list_size"] = long.groupby("query")["db"].transform("size")
    return long


def summarise(rel: pd.DataFrame) -> dict:
    out = dict(n_related_pairs=len(rel), n_queries=rel["query"].nunique(),
               median_rank=float(rel["rank"].median()), median_list_size=float(rel.list_size.median()))
    for k in KS:
        out[f"recall_at_{k}"] = float((rel["rank"] <= k).mean())
        out[f"chance_at_{k}"] = float((np.minimum(k, rel.list_size) / rel.list_size).mean())
    return out


def main():
    pp = pd.read_csv(os.path.join(OUT, "scale_pair_predictions.csv"))
    pp = pp[~pp.same_instrument.astype(bool)].dropna(subset=[VARIANT, "empirical_r"])
    rows, md = [], ["# Rank-based recall of related scales (variant B, cross-instrument)", "",
                    "For each query scale, cross-instrument database scales are ranked by |predicted cosine|. "
                    "Each related pair enters twice (once per direction). Chance = expected share under random "
                    "ordering, min(k, list size) / list size averaged over pairs.", ""]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for ds, d in pp.groupby("dataset"):
        long = ranks_for(d)
        fig, ax = plt.subplots(figsize=(6, 4.2))
        md.append(f"## {ds}")
        md.append("")
        md.append(f"{long['query'].nunique()} query scales, list sizes {int(long.list_size.min())} to {int(long.list_size.max())}.")
        md.append("")
        hdr = ["threshold", "related pairs", "queries", "median rank", "top 5", "top 10", "top 25", "top 50", "top 100"]
        lines = ["| " + " | ".join(hdr) + " |", "|" + "|".join(["---"] * len(hdr)) + "|"]
        for t, color in zip(THRESHOLDS, ["#1b6ca8", "#2a9d5c"]):
            rel = long[long.abs_emp >= t]
            s = summarise(rel)
            rows.append(dict(dataset=ds, threshold=t, **s))
            cells = [f"|r| >= {t:.2f}", str(s["n_related_pairs"]), str(s["n_queries"]), f"{s['median_rank']:.0f}"] + [
                f"{s[f'recall_at_{k}']:.2f} (chance {s[f'chance_at_{k}']:.2f})" for k in KS]
            lines.append("| " + " | ".join(cells) + " |")
            # cumulative recall curve vs rank, and chance
            max_rank = int(long.list_size.max())
            r = np.arange(1, max_rank + 1)
            cum = np.array([(rel["rank"] <= k).mean() for k in r])
            chance = np.array([(np.minimum(k, rel.list_size) / rel.list_size).mean() for k in r])
            ax.plot(r, cum, color=color, label=f"|r| >= {t:.2f} (n = {len(rel)})")
            ax.plot(r, chance, color=color, linestyle="--", linewidth=1, label=f"chance, |r| >= {t:.2f}")
        md.extend(lines)
        md.append("")
        md.append(f"![](rank_recall_{ds}.png)")
        md.append("")
        ax.set_xlabel("rank k in the query's list (by |predicted cosine|)")
        ax.set_ylabel("cumulative recall of related scales")
        ax.set_xscale("log")
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
        ax.set_title(f"Cumulative recall against rank: {ds}")
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"rank_recall_{ds}.png"), dpi=150)
        plt.close(fig)

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "rank_recall.csv"), index=False)
    with open(os.path.join(OUT, "rank_recall.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("\n".join(l for l in md if not l.startswith("![")))


if __name__ == "__main__":
    main()
