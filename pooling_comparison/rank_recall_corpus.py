"""How would rank-based recall look in the real SynthNet corpus?

    python pooling_comparison/rank_recall_corpus.py

The validation datasets have only ~100 candidate scales per query, so rank
recall there says little about a corpus of ~89,000 scales. Here every related
validation pair (query q, database scale d, |empirical r| >= t) is placed into
the real corpus: q is pooled with documented-sign pooling (variant B), scored
against all corpus scale vectors, and the rank d would have in q's corpus
list is 1 + the number of corpus scales whose |cosine| with q exceeds the
predicted |cosine| of (q, d). This assumes d would score in the corpus as it
does here, and treats every corpus scale that outranks d as a distractor,
which is pessimistic: the corpus contains the validation instruments
themselves and many genuinely related scales. A sensitivity variant drops
corpus scales with |cosine| >= .80 to the query (near-duplicates, likely the
same instrument).

Also reported: how many corpus scales exceed the display thresholds per
query, i.e. how long the displayed list would be.

Inputs: results/pooling_comparison/scale_pair_predictions.csv (related pairs and
their predicted cosines), the validation feathers (query vectors) and the
corpus export CORPUS (item embeddings and keying per scale; re-pooled with variant B). Writes
rank_recall_corpus.{csv,md} and rank_recall_corpus_<dataset>.png. Existing
outputs are not modified.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scale_pooling as sp  # noqa: E402
from run_pooling_comparison import DATASETS, Dataset  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "results", "pooling_comparison")
CORPUS = os.environ.get("SYNTHNET_CORPUS", "/Users/rubenarslan/research/items_llm/data/tmp.pkl")
VARIANT = "documented"
THRESHOLDS = [0.40, 0.30]
KS = [5, 10, 25, 50, 100, 250, 1000]
DUP = 0.80


def load_corpus():
    """Re-pool every corpus scale with variant B from its item embeddings and
    keying ('negative' = reverse-keyed; missing keying = positive). The pooled
    vectors shipped in the corpus export are plain means that ignore keying, so
    they cannot be compared with variant-B query vectors."""
    t = pd.read_pickle(CORPUS)[["meta_doi", "is_subscale", "keying", "embeddings"]]
    V = np.empty((len(t), 768), dtype=np.float32)
    n_mismatch = 0
    for k, (E, keys) in enumerate(zip(t.embeddings.values, t.keying.values)):
        E = np.asarray(E, dtype=np.float32)
        if E.ndim == 1:
            E = E[None, :]
        rev = np.array([kk == "negative" for kk in keys], bool) if len(keys) == E.shape[0] else np.zeros(E.shape[0], bool)
        n_mismatch += int(len(keys) != E.shape[0])
        V[k] = sp.pool_documented(E, rev)
    print(f"corpus: {len(t):,} scales re-pooled with variant B; {n_mismatch} rows with keying/embedding length mismatch (treated as all positive)")
    return t.meta_doi.to_numpy(), t.is_subscale.to_numpy(), sp.unit(V)


def main():
    pp = pd.read_csv(os.path.join(OUT, "scale_pair_predictions.csv"))
    pp = pp[~pp.same_instrument.astype(bool)].dropna(subset=[VARIANT, "empirical_r"])
    doi, is_sub, C = load_corpus()
    n_corpus = len(C)
    rows, md = [], [
        "# Rank-based recall extrapolated to the SynthNet corpus (variant B, cross-instrument)", "",
        f"Corpus: {n_corpus:,} scale vectors ({(~is_sub).sum():,} instruments, {is_sub.sum():,} subscales), re-pooled with variant B from item embeddings and keying. "
        "For each related validation pair the rank the database scale would have among all corpus scales is "
        "1 + the number of corpus scales with a higher |cosine| to the query than the pair's predicted |cosine|. "
        f"'excl. near-duplicates' drops corpus scales with |cosine| >= {DUP} to the query. Chance = k / corpus size. "
        "Pessimistic: corpus scales that outrank a related scale may themselves be related to the query.", "",
    ]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for key, tag in DATASETS.items():
        ds = Dataset(key, tag)
        Q = {s: sp.unit(sp.pool_documented(np.stack([ds.emb[i] for i in ds.items[s]]), ds.rev_flags[s])) for s in ds.scale_names}
        names = list(Q)
        S = np.abs(np.stack([Q[s] for s in names]).astype(np.float32) @ C.T)  # queries x corpus
        S_sorted = -np.sort(-S, axis=1)
        qi = {s: k for k, s in enumerate(names)}
        d = pp[pp.dataset == key]
        long = pd.concat([
            d.rename(columns={"scale_a": "query", "scale_b": "db"})[["query", "db", "empirical_r", VARIANT]],
            d.rename(columns={"scale_b": "query", "scale_a": "db"})[["query", "db", "empirical_r", VARIANT]],
        ], ignore_index=True)
        long["abs_pred"] = long[VARIANT].abs()
        long["abs_emp"] = long.empirical_r.abs()
        q_idx = long["query"].map(qi).to_numpy()
        # rank in corpus: number of corpus scales strictly above the pair's cosine
        above = np.array([np.searchsorted(-S_sorted[q], -s, side="left") for q, s in zip(q_idx, long.abs_pred)])
        dup = np.array([(S[q] >= DUP).sum() for q in q_idx])
        long["corpus_rank"] = above + 1
        long["corpus_rank_excl_dup"] = np.maximum(above - dup, 0) + 1
        # displayed list length per query
        shown = {t: (S >= t).sum(axis=1) for t in THRESHOLDS}
        shown_dup = (S >= DUP).sum(axis=1)

        md.append(f"## {key}")
        md.append("")
        md.append(f"{len(names)} query scales. Corpus scales per query with |cosine| >= .40: median {np.median(shown[0.40]):.0f} "
                  f"(IQR {np.percentile(shown[0.40], 25):.0f} to {np.percentile(shown[0.40], 75):.0f}); >= .30: median {np.median(shown[0.30]):.0f}; "
                  f">= {DUP} (near-duplicates): median {np.median(shown_dup):.0f}, max {shown_dup.max()}.")
        md.append("")
        hdr = ["threshold", "ranking", "related pairs", "median corpus rank"] + [f"top {k}" for k in KS]
        lines = ["| " + " | ".join(hdr) + " |", "|" + "|".join(["---"] * len(hdr)) + "|"]
        fig, ax = plt.subplots(figsize=(6, 4.2))
        for t, color in zip(THRESHOLDS, ["#1b6ca8", "#2a9d5c"]):
            rel = long[long.abs_emp >= t]
            for col, label, ls in [("corpus_rank", "all corpus scales", "-"), ("corpus_rank_excl_dup", "excl. near-duplicates", "-.")]:
                r = rel[col].to_numpy()
                row = dict(dataset=key, threshold=t, ranking=label, n_related_pairs=len(rel),
                           n_corpus=n_corpus, median_corpus_rank=float(np.median(r)),
                           median_shown_at_threshold=float(np.median(shown[t])))
                for k in KS:
                    row[f"recall_at_{k}"] = float((r <= k).mean())
                    row[f"chance_at_{k}"] = k / n_corpus
                rows.append(row)
                lines.append("| " + " | ".join(
                    [f"|r| >= {t:.2f}", label, str(len(rel)), f"{np.median(r):.0f}"] +
                    [f"{row[f'recall_at_{k}']:.2f}" for k in KS]) + " |")
                ks = np.unique(np.concatenate([np.arange(1, 101), np.geomspace(100, n_corpus, 60).astype(int)]))
                ax.plot(ks, [(r <= k).mean() for k in ks], color=color, linestyle=ls,
                        label=f"|r| >= {t:.2f}, {label}")
        ks = np.geomspace(1, n_corpus, 100)
        ax.plot(ks, ks / n_corpus, color="grey", linestyle="--", linewidth=1, label="chance")
        ax.set_xscale("log")
        from matplotlib.ticker import FixedLocator, FixedFormatter, NullFormatter
        ticks = [1, 10, 100, 1000, 10000, n_corpus]
        ax.xaxis.set_major_locator(FixedLocator(ticks))
        ax.xaxis.set_major_formatter(FixedFormatter([f"{x:,}" for x in ticks]))
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylim(0, 1)
        ax.grid(alpha=0.3)
        ax.set_xlabel("rank k in the query's corpus list (by |cosine|)")
        ax.set_ylabel("cumulative recall of related scales")
        ax.set_title(f"Recall against corpus rank: {key}")
        ax.legend(frameon=False, fontsize=7, loc="upper left")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"rank_recall_corpus_{key}.png"), dpi=150)
        plt.close(fig)
        md.extend(lines)
        md.append("")
        md.append(f"Chance at k: {', '.join(f'top {k} = {k / n_corpus:.4f}' for k in KS)}.")
        md.append("")
        md.append(f"![](rank_recall_corpus_{key}.png)")
        md.append("")

    pd.DataFrame(rows).to_csv(os.path.join(OUT, "rank_recall_corpus.csv"), index=False)
    with open(os.path.join(OUT, "rank_recall_corpus.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("\n".join(l for l in md if not l.startswith("![")))


if __name__ == "__main__":
    main()
