#!/usr/bin/env python3
"""
STEP 4 - statistics and figures (Figures 1, 2, 4 and the post hoc analyses), as defined in the Methods.
Reads the CSVs written by 03_snp_metrics_and_concordance.py. Needs numpy, scipy, matplotlib.

Usage
  python 04_stats_and_figures.py --results results/ --out results/step3/

Writes
  stats_summary.txt / stats_summary.csv  every number quoted in the Results text
  Figure1_sequencing_metrics.(png|pdf)   Figure2_SNP_count_correlation.(png|pdf)
  Figure4_platform_specific_SNPs.(png|pdf)

Methods followed
  paired two-sided Wilcoxon signed-rank test (n = 21 pairs), alpha 0.05, p < 0.001 reported as such;
            Table 1 p-values are reported raw AND Benjamini-Hochberg adjusted (family = the 7 Table 1 tests)
  Figure 1: boxes = median and IQR, whiskers = most extreme value within 1.5 x IQR, one point per isolate,
            lines join the same isolate across platforms, panels c-e on a log scale
  Figure 2: Pearson, Spearman, OLS slope and R^2, y = x line
  post hoc (exploratory): correlations repeated within the two Illumina-SNP-count groups,
            ONT ambiguous SNPs vs ONT average depth (Spearman), all isolates and the higher-count group
"""
import argparse, csv, os, sys
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ILL_COL, ONT_COL = "#1b6ca8", "#d1603d"
GROUP_SPLIT = 300          # Illumina high-confidence SNPs: low group 169-174, high group 498-656


def read_csv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def load(results):
    t1 = read_csv(os.path.join(results, "table1_per_isolate_metrics.csv"))
    t4 = [r for r in read_csv(os.path.join(results, "table4_recovery_precision.csv"))
          if r["isolate"] not in ("Pooled", "Mean ± SD")]
    isolates = sorted({r["isolate"] for r in t1})
    metrics = {}
    for col in ("avg_depth", "genome_with_coverage_pct", "no_coverage_bases", "no_coverage_pct",
                "high_conf_snps", "ambiguous_snps_AC1", "quality_snps_QUAL_gt150"):
        metrics[col] = {p: np.array([float(next(r[col] for r in t1 if r["isolate"] == i and r["platform"] == p))
                                     for i in isolates]) for p in ("illumina", "ont")}
    t4d = {r["isolate"]: r for r in t4}
    return isolates, metrics, t4d


def fmt_p(p):
    """'<0.001' or '= 0.019' style, ready to follow the letter p."""
    return "< 0.001" if p < 0.001 else f"= {p:.3f}" if p < 0.1 else f"= {p:.2f}"


def msd(x, nd=2):
    return f"{np.mean(x):.{nd}f} ± {np.std(x, ddof=1):.{nd}f} ({np.min(x):.{nd}f}–{np.max(x):.{nd}f})"


# ---------------------------------------------------------------- Figure 1
def paired_panel(ax, a, b, title, ylabel, log=False, p=None):
    data = [a, b]
    bp = ax.boxplot(data, positions=[1, 2], widths=0.5, whis=1.5, showfliers=False,
                    patch_artist=True, medianprops=dict(color="black", lw=1.5))
    for patch, c in zip(bp["boxes"], (ILL_COL, ONT_COL)):
        patch.set_facecolor(c); patch.set_alpha(0.35); patch.set_edgecolor("black")
    for x, y in zip(a, b):
        ax.plot([1, 2], [x, y], color="grey", lw=0.6, alpha=0.7, zorder=2)
    ax.scatter(np.ones(len(a)), a, s=14, color=ILL_COL, zorder=3)
    ax.scatter(np.full(len(b), 2), b, s=14, color=ONT_COL, zorder=3)
    if log:
        ax.set_yscale("log")
    ax.set_xticks([1, 2]); ax.set_xticklabels(["Illumina", "ONT"])
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left", fontsize=10)
    if p is not None:
        ax.text(0.5, 0.98, f"p {fmt_p(p)}", transform=ax.transAxes,
                ha="center", va="top", fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)


def figure1(m, pv, out):
    fig, axs = plt.subplots(2, 3, figsize=(11, 7))
    spec = [("genome_with_coverage_pct", "(a)", "Genome with coverage (%)", False),
            ("avg_depth", "(b)", "Average depth (×)", False),
            ("no_coverage_bases", "(c)", "Bases with zero coverage (n)", True),
            ("no_coverage_pct", "(d)", "Reference with zero coverage (%)", True),
            ("ambiguous_snps_AC1", "(e)", "Ambiguous SNPs (n)", True),
            ("quality_snps_QUAL_gt150", "(f)", "Quality SNPs (n)", False)]
    for ax, (k, lab, yl, log) in zip(axs.ravel(), spec):
        paired_panel(ax, m[k]["illumina"], m[k]["ont"], lab, yl, log, pv[k])
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}/Figure1_sequencing_metrics.{ext}", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- Figure 2
def figure2(x, y, res, out):
    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    ax.scatter(x, y, s=28, color=ONT_COL, edgecolor="black", lw=0.5, zorder=3)
    lim = max(x.max(), y.max()) * 1.05
    ax.plot([0, lim], [0, lim], "--", color="grey", lw=1, label="y = x")
    xs = np.array([0, lim])
    ax.plot(xs, res["intercept"] + res["slope"] * xs, "-", color="black", lw=1.2,
            label=f"least-squares (slope {res['slope']:.2f})")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Illumina high-confidence SNPs (n)"); ax.set_ylabel("ONT high-confidence SNPs (n)")
    ax.text(0.04, 0.96, f"Pearson r = {res['pearson_r']:.3f}\nSpearman ρ = {res['spearman_rho']:.3f}\n"
            f"R² = {res['r2']:.3f}", transform=ax.transAxes, va="top", fontsize=9)
    ax.legend(loc="lower right", frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}/Figure2_SNP_count_correlation.{ext}", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- Figure 4
def figure4(isolates, t4d, out):
    def col(name):
        return np.array([int(t4d[i][name]) for i in isolates])
    panels = [("(a)", "Illumina-only SNPs (n = %d), by ONT call", "ill_only_ONT_ref", "ill_only_ONT_ambig", "ill_only_ONT_uncalled"),
              ("(b)", "ONT-only SNPs (n = %d), by Illumina call", "ont_only_ILL_ref", "ont_only_ILL_ambig", "ont_only_ILL_uncalled")]
    colors = ("#6c8ebf", "#e0a030", "#555555")
    labels = ("Reference", "Ambiguous (AC = 1)", "Uncalled (zero coverage)")
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    x = np.arange(len(isolates))
    for ax, (lab, title, *cols) in zip(axs, panels):
        vals = [col(c) for c in cols]
        bottom = np.zeros(len(isolates))
        for v, c, l in zip(vals, colors, labels):
            ax.bar(x, v, bottom=bottom, color=c, label=l, edgecolor="white", lw=0.4)
            bottom += v
        ax.set_title(f"{lab} " + title % int(sum(v.sum() for v in vals)), loc="left", fontsize=10)
        ax.set_ylabel("SNPs per isolate (n)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper left")
    axs[1].set_xticks(x); axs[1].set_xticklabels(isolates, rotation=60, ha="right", fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}/Figure4_platform_specific_SNPs.{ext}", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="results/step3")
    ap.add_argument("--split", type=int, default=GROUP_SPLIT,
                    help="Illumina SNP count separating the low and high groups")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    isolates, m, t4d = load(a.results)
    n = len(isolates)
    lines, rows = [], []

    def say(s=""):
        lines.append(s); print(s)

    say(f"n isolates = {n}")
    # ---- Table 1 / Figure 1: paired Wilcoxon, descriptive stats, direction counts
    names = {"genome_with_coverage_pct": "Genome with coverage (%)", "avg_depth": "Average depth (x)",
             "no_coverage_bases": "No-coverage bases (n)", "no_coverage_pct": "No-coverage (% of reference)",
             "high_conf_snps": "High-confidence SNPs (n)", "ambiguous_snps_AC1": "Ambiguous SNPs (n)",
             "quality_snps_QUAL_gt150": "Quality SNPs (n)"}
    pv = {}
    say("Wilcoxon: scipy default zero_method='wilcox' (isolates with identical values on both platforms are dropped); "
        "exact p-value when n<=50 with no ties/zeros, otherwise normal approximation.")
    say("\n== Paired Wilcoxon signed-rank (two-sided), Illumina vs ONT ==")
    for k, label in names.items():
        i_, o_ = m[k]["illumina"], m[k]["ont"]
        res = stats.wilcoxon(i_, o_, alternative="two-sided")
        pv[k] = res.pvalue
        nd = 0 if "bases" in k or k.endswith("snps") or "SNPs" in label else 2
        nd = 1 if k == "ambiguous_snps_AC1" else nd
        higher, lower, ties = int((o_ > i_).sum()), int((o_ < i_).sum()), int((o_ == i_).sum())
        say(f"{label:32s} Illumina {msd(i_, nd)} | ONT {msd(o_, nd)} | median {np.median(i_):.2f} vs "
            f"{np.median(o_):.2f} | ONT higher in {higher}, lower in {lower}, equal in {ties} | "
            f"W = {res.statistic:.1f}, p {fmt_p(res.pvalue)} (raw {res.pvalue:.4g})")
        rows.append([label, msd(i_, nd), msd(o_, nd), f"{np.median(i_):.2f}", f"{np.median(o_):.2f}",
                     higher, lower, ties, f"{res.statistic:.1f}", f"{res.pvalue:.4g}"])
    # Benjamini-Hochberg over the Table 1 family
    ks = list(names)
    padj = stats.false_discovery_control([pv[k] for k in ks], method="bh")
    say("\n== Benjamini-Hochberg adjustment (family = the %d Table 1 tests) ==" % len(ks))
    for k, pa in zip(ks, padj):
        say(f"{names[k]:32s} raw p {fmt_p(pv[k])} (raw {pv[k]:.4g}) | BH-adjusted p {fmt_p(pa)} (adj {pa:.4g})")
    sig_raw = [k for k in ks if pv[k] < 0.05]
    lost = [names[k] for k, pa in zip(ks, padj) if pv[k] < 0.05 and pa >= 0.05]
    say(f"raw p < 0.05: {len(sig_raw)} of {len(ks)}; still < 0.05 after BH: {len(sig_raw) - len(lost)}"
        + (f"; LOST significance: {lost}  --> edit the Methods sentence accordingly" if lost else
           "; every p < 0.05 remained < 0.05 (supports the Methods sentence)"))
    for k in ("no_coverage_bases", "no_coverage_pct", "ambiguous_snps_AC1"):
        nz = int((m[k]["illumina"] <= 0).sum() + (m[k]["ont"] <= 0).sum())
        if nz:
            say(f"WARNING: {nz} zero values in '{k}' (log-scale panel cannot plot zeros; state this in the legend)")
    # extra figures quoted in the text
    cov_i, cov_o = m["genome_with_coverage_pct"]["illumina"], m["genome_with_coverage_pct"]["ont"]
    nc_o = m["no_coverage_bases"]["ont"]
    hi_cov = cov_o >= 99.9
    def rng(v, f="{:.0f}"):
        return f"{f.format(v.min())}-{f.format(v.max())}" if len(v) else "none"
    say(f"\nONT isolates with >=99.9% coverage: {int(hi_cov.sum())} (uncovered {rng(nc_o[hi_cov])}); "
        f"others: {int((~hi_cov).sum())} (coverage {rng(cov_o[~hi_cov], '{:.2f}')}%, uncovered {rng(nc_o[~hi_cov])})")
    d_o = m["avg_depth"]["ont"]
    say(f"ONT depth range {d_o.min():.1f}-{d_o.max():.1f}; isolates with ONT depth <13.8x (post hoc threshold, exploratory): {int((d_o < 13.8).sum())}; "
        f"Illumina depth range {m['avg_depth']['illumina'].min():.1f}-{m['avg_depth']['illumina'].max():.1f}")
    say(f"No-coverage %: Illumina {m['no_coverage_pct']['illumina'].min():.2f}-{m['no_coverage_pct']['illumina'].max():.2f}, "
        f"ONT {m['no_coverage_pct']['ont'].min():.3f}-{m['no_coverage_pct']['ont'].max():.2f}")
    amb_o = m["ambiguous_snps_AC1"]["ont"]
    say(f"ONT ambiguous max {int(amb_o.max())} at depth {d_o[amb_o.argmax()]:.1f}x (isolate {isolates[int(amb_o.argmax())]}); "
        f"isolate with lowest ONT depth: {isolates[int(d_o.argmin())]} ({d_o.min():.1f}x)")

    # ---- Figure 2: SNP count correlation
    x, y = m["high_conf_snps"]["illumina"], m["high_conf_snps"]["ont"]
    lr = stats.linregress(x, y)
    res2 = dict(pearson_r=stats.pearsonr(x, y)[0], spearman_rho=stats.spearmanr(x, y)[0], r2=lr.rvalue ** 2,
                slope=lr.slope, intercept=lr.intercept)
    say("\n== SNP count correlation (all isolates) ==")
    say(f"Pearson r = {res2['pearson_r']:.3f}, Spearman rho = {res2['spearman_rho']:.3f}, R2 = {res2['r2']:.3f}, "
        f"slope = {res2['slope']:.2f}, intercept = {res2['intercept']:.1f}; ONT > Illumina in {int((y > x).sum())}/{n}")
    say("\n== Post hoc (exploratory, not corrected for multiple comparisons) ==")
    hi = x >= a.split
    say(f"groups by Illumina SNP count (split at {a.split}): low n={int((~hi).sum())} "
        f"({int(x[~hi].min())}-{int(x[~hi].max())}), high n={int(hi.sum())} ({int(x[hi].min())}-{int(x[hi].max())})")
    for name, mask in (("low group", ~hi), ("high group", hi)):
        r, p = stats.pearsonr(x[mask], y[mask]); rho, ps = stats.spearmanr(x[mask], y[mask])
        say(f"{name:10s} Pearson r = {r:.2f} (p {fmt_p(p)}), Spearman rho = {rho:.2f} (p {fmt_p(ps)}), n = {int(mask.sum())}")
    rho, p = stats.spearmanr(amb_o, d_o)
    say(f"ONT ambiguous vs ONT depth, all {n}: Spearman rho = {rho:.2f}, p {fmt_p(p)} (raw {p:.4g})")
    rho, p = stats.spearmanr(amb_o[hi], d_o[hi])
    say(f"ONT ambiguous vs ONT depth, high group n={int(hi.sum())}: Spearman rho = {rho:.2f}, p {fmt_p(p)} (raw {p:.4g})")

    # ---- Table 4 summaries quoted in the text (computed from table4 CSV)
    say("\n== Recovery / precision (from table4) ==")
    rec = np.array([float(t4d[i]["recovery_pct"]) for i in isolates])
    prec = np.array([float(t4d[i]["ont_precision_pct"]) for i in isolates])
    say(f"recovery {rec.mean():.1f} ± {rec.std(ddof=1):.1f} (range {rec.min():.1f}-{rec.max():.1f}); "
        f"precision {prec.mean():.1f} ± {prec.std(ddof=1):.1f} (range {prec.min():.1f}-{prec.max():.1f}; "
        f"lowest {isolates[int(prec.argmin())]}, highest {isolates[int(prec.argmax())]})")
    ont_only = np.array([int(t4d[i]["ont_only"]) for i in isolates])
    ill_only = np.array([int(t4d[i]["illumina_only"]) for i in isolates])
    say(f"Illumina-only per isolate {ill_only.min()}-{ill_only.max()} (total {ill_only.sum()}); "
        f"ONT-only total {ont_only.sum()}; low group {ont_only[~hi].min()}-{ont_only[~hi].max()}, "
        f"high group {ont_only[hi].min()}-{ont_only[hi].max()}")

    # ---- figures
    figure1(m, pv, a.out)
    figure2(x, y, res2, a.out)
    figure4(isolates, t4d, a.out)
    with open(f"{a.out}/stats_summary.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(f"{a.out}/stats_summary.csv", "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "illumina_mean_sd_range", "ont_mean_sd_range", "illumina_median", "ont_median",
                    "ont_higher_n", "ont_lower_n", "equal_n", "wilcoxon_W", "p_raw"])
        w.writerows(rows)
    print(f"\nwrote figures and stats_summary.* to {a.out}")


if __name__ == "__main__":
    main()
