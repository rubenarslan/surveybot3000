"""Threshold-based evaluation of documented-sign pooling (variant B) on
cross-instrument scale pairs, from the existing scale_pair_predictions.csv.

    python pooling_comparison/threshold_evaluation.py

Nothing is recomputed: reads results/pooling_comparison/scale_pair_predictions.csv
and writes threshold_metrics.csv, calibration.csv, calibration_<dataset>.png and
threshold_summary.md next to it. Existing outputs are not modified.

1. Precision and recall at display thresholds t in {.30, .40, .50}, with counts
   and bootstrap 95% intervals over scale pairs, plus the base rate.
2. Calibration: pairs binned by |predicted| (0-.1, ..., .6-.7, .7+), with n,
   mean and median |empirical| and the share with |empirical| >= .40.
3. Sign agreement among pairs with |predicted| >= t (t = .40 is the one that
   matters; the others are reported for context).
4. Precision at 5 by absolute cosine, cross-instrument, restricted to query
   scales with at least one cross-instrument scale at |empirical| >= .40.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results", "pooling_comparison")
VARIANT = "documented"
THRESHOLDS = [0.30, 0.40, 0.50]
BINS = [0, .1, .2, .3, .4, .5, .6, .7, np.inf]
BIN_LABELS = ["0-.1", ".1-.2", ".2-.3", ".3-.4", ".4-.5", ".5-.6", ".6-.7", ".7+"]
RELEVANT = 0.40
N_BOOT = 2000
SEED = 20260925


def precision_recall(pred_abs, emp_abs, t):
    flagged = pred_abs >= t
    related = emp_abs >= t
    both = flagged & related
    prec = both.sum() / flagged.sum() if flagged.sum() else np.nan
    rec = both.sum() / related.sum() if related.sum() else np.nan
    return prec, rec, int(flagged.sum()), int(related.sum()), int(both.sum())


def bootstrap(pred_abs, emp_abs, t, rng):
    n = len(pred_abs)
    ps, rs = [], []
    idx = rng.integers(0, n, size=(N_BOOT, n))
    for i in idx:
        p, r, *_ = precision_recall(pred_abs[i], emp_abs[i], t)
        ps.append(p)
        rs.append(r)
    ps, rs = np.array(ps, float), np.array(rs, float)
    return (np.nanpercentile(ps, 2.5), np.nanpercentile(ps, 97.5),
            np.nanpercentile(rs, 2.5), np.nanpercentile(rs, 97.5))


def precision_at_5(d):
    long = pd.concat([
        d.rename(columns={"scale_a": "query", "scale_b": "db"})[["query", "db", "empirical_r", VARIANT]],
        d.rename(columns={"scale_b": "query", "scale_a": "db"})[["query", "db", "empirical_r", VARIANT]],
    ])
    vals = []
    for _, g in long.groupby("query"):
        rel = np.abs(g.empirical_r.to_numpy()) >= RELEVANT
        if not rel.any():
            continue
        order = np.argsort(-np.abs(g[VARIANT].to_numpy()), kind="stable")
        vals.append(rel[order][:5].mean())
    return float(np.mean(vals)), len(vals)


def md_table(df, cols, fmt):
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmt.get(c, "{}").format(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def main():
    pp = pd.read_csv(os.path.join(OUT, "scale_pair_predictions.csv"))
    pp = pp[~pp.same_instrument.astype(bool)].dropna(subset=[VARIANT, "empirical_r"])
    rng = np.random.default_rng(SEED)
    thr_rows, cal_rows, md = [], [], []
    md.append("# Threshold-based evaluation of documented-sign pooling (variant B)")
    md.append("")
    md.append("Cross-instrument scale pairs only (query and database scale from different instruments), "
              "from `scale_pair_predictions.csv`. Bootstrap intervals: 2000 resamples of scale pairs.")
    md.append("")
    for ds, d in pp.groupby("dataset"):
        pred, emp = d[VARIANT].to_numpy(), d.empirical_r.to_numpy()
        pa, ea = np.abs(pred), np.abs(emp)
        p5, nq = precision_at_5(d)
        md.append(f"## {ds}")
        md.append("")
        md.append(f"{len(d)} cross-instrument pairs among {len(set(d.scale_a) | set(d.scale_b))} scales.")
        md.append("")
        for t in THRESHOLDS:
            prec, rec, nf, nr, nb = precision_recall(pa, ea, t)
            plo, phi, rlo, rhi = bootstrap(pa, ea, t, rng)
            flagged = pa >= t
            sign_agree = float(np.mean(np.sign(pred[flagged]) == np.sign(emp[flagged]))) if flagged.any() else np.nan
            thr_rows.append(dict(
                dataset=ds, threshold=t, n_pairs=len(d), base_rate=float((ea >= t).mean()),
                n_flagged=nf, n_related=nr, n_both=nb,
                precision=prec, precision_lo=plo, precision_hi=phi,
                recall=rec, recall_lo=rlo, recall_hi=rhi,
                sign_agreement_flagged=sign_agree,
                p_at_5_abs=p5, n_queries_with_relevant=nq,
            ))
        # calibration
        b = np.asarray(pd.cut(pa, BINS, labels=BIN_LABELS, right=False, include_lowest=True).astype(str))
        for lab in BIN_LABELS:
            m = b == lab
            cal_rows.append(dict(
                dataset=ds, bin=lab, n=int(m.sum()),
                mean_abs_empirical=float(ea[m].mean()) if m.any() else np.nan,
                median_abs_empirical=float(np.median(ea[m])) if m.any() else np.nan,
                share_empirical_ge_40=float((ea[m] >= RELEVANT).mean()) if m.any() else np.nan,
                mean_abs_predicted=float(pa[m].mean()) if m.any() else np.nan,
            ))
        # plot
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cal = pd.DataFrame([r for r in cal_rows if r["dataset"] == ds and r["n"] > 0])
        fig, ax = plt.subplots(figsize=(5.5, 5))
        ax.scatter(pa, ea, s=6, alpha=0.25, color="#1b6ca8", linewidths=0)
        ax.plot(cal.mean_abs_predicted, cal.mean_abs_empirical, marker="o", color="#d1495b", label="bin mean")
        ax.plot([0, 1], [0, 1], color="grey", linestyle="--", linewidth=1, label="identity")
        ax.axhline(RELEVANT, color="#2a9d5c", linewidth=0.8, linestyle=":")
        ax.axvline(RELEVANT, color="#2a9d5c", linewidth=0.8, linestyle=":")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.set_xlabel("|predicted| (documented-sign cosine)")
        ax.set_ylabel("|empirical r|")
        ax.set_title(f"Calibration, cross-instrument pairs: {ds}")
        ax.legend(frameon=False, loc="upper left")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"calibration_{ds}.png"), dpi=150)
        plt.close(fig)

        tr = pd.DataFrame([r for r in thr_rows if r["dataset"] == ds])
        tr["precision_ci"] = tr.apply(lambda r: f"{r.precision:.2f} [{r.precision_lo:.2f}, {r.precision_hi:.2f}]", axis=1)
        tr["recall_ci"] = tr.apply(lambda r: f"{r.recall:.2f} [{r.recall_lo:.2f}, {r.recall_hi:.2f}]", axis=1)
        md.append("### Precision and recall at display threshold t")
        md.append("")
        md.append(md_table(tr, ["threshold", "base_rate", "n_flagged", "n_related", "n_both", "precision_ci", "recall_ci", "sign_agreement_flagged"],
                           {"threshold": "{:.2f}", "base_rate": "{:.3f}", "sign_agreement_flagged": "{:.3f}"}))
        md.append("")
        md.append("### Calibration by |predicted| bin")
        md.append("")
        md.append(md_table(cal, ["bin", "n", "mean_abs_empirical", "median_abs_empirical", "share_empirical_ge_40"],
                           {"mean_abs_empirical": "{:.3f}", "median_abs_empirical": "{:.3f}", "share_empirical_ge_40": "{:.3f}"}))
        md.append("")
        md.append(f"![]({os.path.basename(os.path.join(OUT, f'calibration_{ds}.png'))})")
        md.append("")
        md.append(f"Precision at 5 by |cosine|, cross-instrument, queries with at least one related scale: "
                  f"{p5:.3f} over {nq} query scales.")
        md.append("")
    pd.DataFrame(thr_rows).to_csv(os.path.join(OUT, "threshold_metrics.csv"), index=False)
    pd.DataFrame(cal_rows).to_csv(os.path.join(OUT, "calibration.csv"), index=False)
    with open(os.path.join(OUT, "threshold_summary.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
