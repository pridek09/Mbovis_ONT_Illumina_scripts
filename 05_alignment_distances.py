#!/usr/bin/env python3
"""
STEP 5 - pairwise SNP distances from the vSNP3 step 2 alignments, three masking conditions, Figure 3.
This is the PRIMARY distance analysis described in the Methods section "Pairwise SNP distance".

Inputs (per platform): the step 2 sorted_table .xlsx and the step 2 SNP alignment .fasta (Mbovis-All group),
plus the duplicate-removed BAMs of step 1 (for the callability filter).

Conditions
  (i)   unfiltered         : the vSNP3 step 2 alignment as written (no additional masking)
  (ii)  PE_PPE_only        : sites whose vSNP3 annotation is PE / PPE / PE_PGRS removed
  (iii) PE_PPE_plus_callability : (ii) plus the read-level callability filter below

Callability filter - a site is removed from BOTH platforms' alignments if, on either platform:
  R2a  the median depth over isolates that have reads at the site is < --min-depth (default 5)
  R2b  no isolate has any read at the site (zero coverage in every isolate)
  R3   more than --max-lowmq (default 0.25) of the pooled reads at the site have MQ < 20
  R4   the median depth in that platform's isolates carrying the alternative base is < --min-depth
Ambiguity codes alone are NOT a criterion (an ambiguous call can be a zero-coverage isolate).

Distance: number of alignment positions at which both sequences carry A/C/G/T and the bases differ
(identical to snp-dists). The root sequence is excluded from pairs. The unfiltered Python matrix is checked
against the snp-dists TSV next to each FASTA when present.

Outputs (in --out)
  comparison_summary.tsv, distance_summary.txt, removed_sites_report.tsv, pairwise_comparison_<condition>.csv,
  <Platform>_<condition>_snp_distances.tsv, <Platform>_<condition>.fasta, Figure3_pairwise_SNP_distance.(png|pdf)
  depth_ill.tsv / depth_ont.tsv are cached; the cache is rebuilt automatically when the BAM list or sites change.

Condition (iv) - OPTIONAL, added with --curated-xlsx: the curated matrices (Final_curated_SNPs_ONT_Illumina_v6.xlsx,
sheets Cleaned_ONT and Cleaned_Illumina). They already contain every filtering decision (PE/PPE, <=10 bp clusters,
IGV review), so NO extra masking is applied to them (the PE/PPE and callability masks would remove the clusters you
validated and kept). Three curated conditions are written:
  curated                  each platform's own cleaned matrix (ONT 833 sites, Illumina 824)
  curated_shared_sites     the sites present in both matrices
  curated_complete_case    shared sites that are A/C/G/T in all 21 isolates on both platforms
Rows are matched by isolate NAME (ONT names with two ERR accessions are reduced to the first), checked against
--samples; N, IUPAC codes and '-' are never counted as differences (pairwise deletion, as in snp-dists).

Run (all conditions)
  python3 05_alignment_distances.py --ill-xlsx ... --ill-fasta ... --ont-xlsx ... --ont-fasta ... \
      --curated-xlsx Final_curated_SNPs_ONT_Illumina_v6.xlsx --samples samples.csv --out OUT
Run (curated only, no BAMs/FASTA needed)
  python3 05_alignment_distances.py --curated-only --curated-xlsx Final_curated_SNPs_ONT_Illumina_v6.xlsx --samples samples.csv --out OUT
"""
import argparse, collections, concurrent.futures as cf, glob, hashlib, itertools, os, re, subprocess, sys
import numpy as np, openpyxl, pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = "/home/pridek/Illumina_ONT/analyses/vsnp3_supp"
TAG = {"ill": "Illumina", "ont": "ONT"}
# PE/PPE/PE_PGRS annotation (word boundaries on every alternative so "type family" etc. do not match)
PAT = re.compile(r"PE[_-]?PGRS|\bPPE|\bPE\d|\bPE_|\bPE-|\bppe family|\bpe family", re.I)


def norm(x):
    """Isolate name: drop run prefix, _nodup and platform/depth suffixes."""
    x = re.sub(r"^\d{2}-\d{6}-", "", str(x).strip())
    x = re.sub(r"(_nodup|_15x|_nanopore|_ont)+$", "", x, flags=re.I)
    m = re.match(r"^(ERR\d+)(?:ERR\d+)+$", x)      # curated ONT sheet: two glued accessions -> first one
    return m.group(1) if m else x


def load_table(f):
    ws = openpyxl.load_workbook(f).active
    rows = [[c.value for c in r] for r in ws.iter_rows()]
    pos = [int(h.split(":")[1]) for h in rows[0][1:]]
    body = rows[1:-2]                       # sample rows (last two rows: summary row and annotation row)
    gt = pd.DataFrame([r[1:] for r in body], index=[norm(r[0]) for r in body], columns=pos)
    return gt, [x or "" for x in rows[-1][1:]], [r[0] for r in body]


def read_fasta(f):
    d, k = {}, None
    for l in open(f):
        l = l.strip()
        if l.startswith(">"):
            k = l[1:].split()[0]; d[k] = []
        elif k:
            d[k].append(l)
    return {k: "".join(v) for k, v in d.items()}


def write_fasta(d, f):
    with open(f, "w") as o:
        for k, s in d.items():
            o.write(f">{k}\n")
            for i in range(0, len(s), 80):
                o.write(s[i:i + 80] + "\n")


def dist(seqs):
    n = list(seqs); M = np.array([list(seqs[x]) for x in n]); ok = np.isin(M, list("ACGT"))
    D = np.zeros((len(n),) * 2, int)
    for i in range(len(n)):
        D[i] = ((M[i] != M) & ok[i] & ok).sum(1)
    return pd.DataFrame(D, n, n)


def read_snpdists(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = [norm(x) for x in df.index]; df.columns = [norm(x) for x in df.columns]
    return df


def summarise(Dill, Dont, label, a):
    """Dill/Dont = (distance DataFrame, number of sites). Writes pairwise_comparison_<label>.csv."""
    ids = sorted((set(Dill[0].index) & set(Dont[0].index)) - {"root"})
    if len(ids) != a.expect:
        sys.exit(f"expected {a.expect} matched isolates, found {len(ids)}")
    c = pd.DataFrame([(x, y, Dill[0].loc[x, y], Dont[0].loc[x, y]) for x, y in itertools.combinations(ids, 2)],
                     columns=["a", "b", "illumina", "ont"])
    c["diff"] = c.ont - c.illumina
    c.to_csv(f"{a.out}/pairwise_comparison_{label}.csv", index=False)
    ad = c["diff"].abs(); cl = c[c.illumina <= a.close]
    s = dict(condition=label, sites_Illumina=Dill[1], sites_ONT=Dont[1], pairs=len(c),
             exact_pct=round(100 * (ad == 0).mean(), 1), le1_pct=round(100 * (ad <= 1).mean(), 1),
             le2_pct=round(100 * (ad <= 2).mean(), 1), le5_pct=round(100 * (ad <= 5).mean(), 1),
             median_diff=float(c["diff"].median()), IQR=f"{c['diff'].quantile(.25):.0f} to {c['diff'].quantile(.75):.0f}",
             range=f"{c['diff'].min()}..{c['diff'].max()}", mean_diff=round(c["diff"].mean(), 2),
             spearman=round(c.illumina.corr(c.ont, method="spearman"), 4),
             close_pairs=len(cl), close_mean_Illumina=round(cl.illumina.mean(), 2), close_mean_ONT=round(cl.ont.mean(), 2))
    return s, c


def load_curated(xlsx, sheet):
    ws = openpyxl.load_workbook(xlsx, data_only=True)[sheet]
    rows = [list(r) for r in ws.iter_rows(values_only=True)]
    if rows[-2][0] != "MQ" or rows[-1][0] != "annotation":
        sys.exit(f"{sheet}: last two rows are not MQ / annotation")
    pos = [int(str(h).split(":")[1]) for h in rows[0][1:]]
    body = rows[1:-2]
    names = ["root" if r[0] == "root" else norm(r[0]) for r in body]
    if len(set(names)) != len(names) or len(set(pos)) != len(pos):
        sys.exit(f"{sheet}: duplicate isolate names or positions after harmonising")
    return pd.DataFrame([[(c or "N").strip().upper() for c in r[1:]] for r in body], index=names, columns=pos)


def curated_conditions(a):
    """Returns {label: ((Dill, n_sites), (Dont, n_sites))} for the curated matrices."""
    import csv
    ont = load_curated(a.curated_xlsx, "Cleaned_ONT"); ill = load_curated(a.curated_xlsx, "Cleaned_Illumina")
    ids = sorted((set(ont.index) & set(ill.index)) - {"root"})
    for lab, g in (("ONT", ont), ("Illumina", ill)):
        extra = sorted(set(g.index) - {"root"} - set(ids))
        if extra or len(ids) != a.expect:
            sys.exit(f"curated {lab}: isolate names do not match the other platform / --expect {a.expect}. Unmatched: {extra}")
    if a.samples:
        exp = sorted({r["isolate"] for r in csv.DictReader(open(a.samples))})
        if exp != ids:
            sys.exit(f"curated isolates differ from {a.samples}: {sorted(set(ids) ^ set(exp))}")
    shared = sorted(set(ont.columns) & set(ill.columns))
    print(f"\nCurated: ONT {ont.shape[1]} sites | Illumina {ill.shape[1]} | shared {len(shared)} | "
          f"ONT-only {ont.shape[1] - len(shared)} | Illumina-only {ill.shape[1] - len(shared)}")
    O, I = ont.loc[ids, shared].values, ill.loc[ids, shared].values
    bo, bi = np.isin(O, list("ACGT")), np.isin(I, list("ACGT")); both = bo & bi; disc = both & (O != I)
    print(f"Curated cells at shared sites: {O.size} | both unambiguous {both.sum()} | discordant {disc.sum()} | "
          f"ambiguous ONT {(~bo).sum()}, Illumina {(~bi).sum()} | identical incl. ambiguity codes {100 * (O == I).mean():.2f}%")
    for r, c in zip(*np.where(disc)):
        print("  DISCORDANT", ids[r], shared[c], "ONT", O[r, c], "Illumina", I[r, c])
    cc = [p for p, ok in zip(shared, both.all(axis=0)) if ok]
    sets = {"curated": (list(ill.columns), list(ont.columns)), "curated_shared_sites": (shared, shared),
            "curated_complete_case": (cc, cc)}
    out = {}
    for label, (pi, po) in sets.items():
        Dp = []
        for gt, pp, tag in ((ill, pi, "ill"), (ont, po, "ont")):
            d = dist({n: "".join(gt.loc[n, pp]) for n in ids})
            d.to_csv(f"{a.out}/{TAG[tag]}_{label}_snp_distances.tsv", sep="\t")
            Dp.append((d, len(pp)))
        out[label] = tuple(Dp)
    return out


def figure3(tabs, res, titles, a):
    labels = list(tabs); n = len(labels); nc = min(n, 3); nr = -(-n // nc)
    fig, axs = plt.subplots(nr, nc, figsize=(4.7 * nc, 4.8 * nr), squeeze=False)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    for ax, label in zip(axs.ravel(), labels):
        c = tabs[label]; x, y = c.illumina.values.astype(float), c.ont.values.astype(float)
        lim = max(x.max(), y.max()) * 1.05 + 1
        ax.plot([0, lim], [0, lim], "--", color="grey", lw=1, label="y = x")
        ax.scatter(x, y, s=18, color="#1b6ca8", edgecolor="white", lw=0.4, zorder=3)
        s = res[res.condition == label].iloc[0]
        ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
        ax.set_xlabel("Illumina pairwise SNP distance"); ax.set_ylabel("ONT pairwise SNP distance")
        ax.set_title(titles.get(label, label), loc="left", fontsize=10)
        ax.text(0.97, 0.03, f"{len(c)} pairs\nSpearman ρ = {s.spearman:.3f}\nexact {s.exact_pct}% | ≤2 {s.le2_pct}% | ≤5 {s.le5_pct}%",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}/Figure3_pairwise_SNP_distance.{ext}", dpi=300)
    plt.close(fig)


def run_raw(a):
    """Conditions (i)-(iii) from the raw vSNP3 step 2 tables/alignments (logic unchanged)."""
    T = {k: load_table(getattr(a, f"{k}_xlsx")) for k in ("ill", "ont")}
    for k in T:
        print(f"{TAG[k]}: {T[k][0].shape[1]} variable sites in the step 2 table, {len(T[k][0]) - 1} samples + root")
    allpos = sorted(set(T["ill"][0].columns) | set(T["ont"][0].columns))

    # ---------------- Step 1: depth at every alignment site (cached, cache validated)
    bed = f"{a.out}/sites.bed"
    open(bed, "w").write("".join(f"{a.chr}\t{p - 1}\t{p}\n" for p in allpos))

    def depth_one(args):
        bam, q = args
        r = subprocess.run(["samtools", "depth", "-a", "-b", bed, "-q", "0", "-Q", str(q), "-d", "0", bam],
                           capture_output=True, text=True)
        if r.returncode:
            sys.exit(f"samtools depth failed on {bam}: {r.stderr[:300]}")
        return [(int(l.split("\t")[1]), int(l.split("\t")[2])) for l in r.stdout.splitlines() if l]

    for k in ("ill", "ont"):
        cache, key = f"{a.out}/depth_{k}.tsv", f"{a.out}/depth_{k}.key"
        bams = sorted(glob.glob(getattr(a, f"{k}_bams")))
        if not bams:
            sys.exit(f"no BAMs for {TAG[k]}: {getattr(a, f'{k}_bams')}")
        sig = hashlib.md5(("|".join(f"{b}:{os.path.getmtime(b)}:{os.path.getsize(b)}" for b in bams)
                           + f"|{len(allpos)}|{allpos[0]}|{allpos[-1]}").encode()).hexdigest()
        if os.path.exists(cache) and os.path.exists(key) and open(key).read() == sig:
            print("using validated depth cache", cache); continue
        print(f"{TAG[k]}: {len(bams)} BAMs, {len(allpos)} sites ...")
        with cf.ThreadPoolExecutor(a.threads) as ex:
            res = list(ex.map(depth_one, [(b, q) for b in bams for q in (0, 20)]))
        rows = []
        for i, b in enumerate(bams):
            d0, d20 = dict(res[2 * i]), dict(res[2 * i + 1])
            rows += [(norm(os.path.basename(b)[:-4]), p, d0.get(p, 0), d20.get(p, 0)) for p in allpos]
        pd.DataFrame(rows, columns=["sample", "pos", "dall", "dq20"]).to_csv(cache, sep="\t", index=False)
        open(key, "w").write(sig)
    DP = {k: pd.read_csv(f"{a.out}/depth_{k}.tsv", sep="\t") for k in ("ill", "ont")}

    # ---------------- Step 2: flag sites
    flag = collections.defaultdict(set)
    for k in ("ill", "ont"):
        gt, ann, _ = T[k]
        for p, t in zip(gt.columns, ann):
            if PAT.search(str(t)):
                flag[p].add("R1_PE_PPE")
        d = DP[k]; cov = d[d.dall > 0]
        med = cov.groupby("pos").dall.median()
        tot = d.groupby("pos")[["dall", "dq20"]].sum()
        lowmq = 1 - tot.dq20 / tot.dall.replace(0, np.nan)
        for p in allpos:
            if p in med.index:
                if med[p] < a.min_depth:
                    flag[p].add(f"R2a_lowdepth_{TAG[k]}")
            else:
                flag[p].add(f"R2b_nocoverage_{TAG[k]}")
            if lowmq.get(p, 0) > a.max_lowmq:
                flag[p].add(f"R3_multimap_{TAG[k]}")
        root = gt.loc["root"]; dd = {(s, p): v for s, p, v in zip(d["sample"], d["pos"], d["dall"])}
        for p in gt.columns:
            carriers = [s for s in gt.index if s != "root" and gt.loc[s, p] in list("ACGT") and gt.loc[s, p] != root[p]]
            v = [dd[(s, p)] for s in carriers if (s, p) in dd]
            if v and np.median(v) < a.min_depth:
                flag[p].add(f"R4_weak_support_{TAG[k]}")
    rep = pd.DataFrame([(p, ";".join(sorted(r))) for p, r in sorted(flag.items())], columns=["pos", "reasons"])
    rep["in_Illumina"] = rep.pos.isin(T["ill"][0].columns); rep["in_ONT"] = rep.pos.isin(T["ont"][0].columns)
    rep.to_csv(f"{a.out}/removed_sites_report.tsv", sep="\t", index=False)
    print(f"\nsites flagged: {len(rep)} (in Illumina table {rep.in_Illumina.sum()}, in ONT table {rep.in_ONT.sum()})")
    for r in sorted({x.split('_')[0] for rs in flag.values() for x in rs}):
        print(f"  {r}: {sum(any(x.startswith(r) for x in rs) for rs in flag.values())} sites")

    # ---------------- Step 3: masking, distances, comparison
    def build(maskpos, label):
        D = {}
        for k in ("ill", "ont"):
            gt, ann, names = T[k]; seqs = read_fasta(getattr(a, f"{k}_fasta"))
            if {len(s) for s in seqs.values()} != {gt.shape[1]}:
                sys.exit(f"{TAG[k]}: FASTA length != number of table sites")
            sm = [n for n in names if n in seqs]
            if len(sm) < 2:
                sys.exit(f"{TAG[k]}: sample names in table and FASTA do not overlap")
            tr = {n: gt.loc[norm(n)].tolist() for n in sm}
            pool = collections.defaultdict(collections.deque)      # table and FASTA column order can differ
            for j in range(gt.shape[1]):
                pool[tuple(seqs[n][j] for n in sm)].append(j)
            drop = set()
            for i, p in enumerate(gt.columns):
                sig = tuple(tr[n][i] for n in sm)
                if not pool[sig]:
                    sys.exit(f"{TAG[k]}: table column {p} not found in FASTA")
                j = pool[sig].popleft()
                if p in maskpos:
                    drop.add(j)
            keep = [j for j in range(gt.shape[1]) if j not in drop]
            new = {n: "".join(s[j] for j in keep) for n, s in seqs.items()}
            d = dist(new)
            d.to_csv(f"{a.out}/{TAG[k]}_{label}_snp_distances.tsv", sep="\t")
            write_fasta(new, f"{a.out}/{TAG[k]}_{label}.fasta")
            d.index = d.columns = [norm(x) for x in d.index]
            D[k] = (d, len(keep))
            if label == "unfiltered":     # cross-check against snp-dists output stored next to the FASTA
                ref = re.sub(r"\.fasta$", "_snp_distances.tsv", getattr(a, f"{k}_fasta"))
                if os.path.exists(ref):
                    sd = read_snpdists(ref); common = [x for x in d.index if x in sd.index]
                    same = bool((sd.loc[common, common].values == d.loc[common, common].values).all())
                    print(f"{TAG[k]}: unfiltered Python distances {'MATCH' if same else 'DO NOT MATCH'} snp-dists ({ref})")
                    if not same:
                        sys.exit("distance calculation disagrees with snp-dists; investigate before using results")
                else:
                    print(f"{TAG[k]}: no snp-dists file found next to FASTA; cross-check skipped")
        return summarise(D["ill"], D["ont"], label, a)

    pe = set(rep[rep.reasons.str.contains("R1")].pos); allflag = set(rep.pos)
    conds = [("unfiltered", set()), ("PE_PPE_only", pe), ("PE_PPE_plus_callability", allflag)]
    summ, tabs = [], {}
    for label, mp in conds:
        s, c = build(mp, label); summ.append(s); tabs[label] = c
    return summ, tabs, conds


def main():
    ap = argparse.ArgumentParser()
    for k in ("ill", "ont"):
        ap.add_argument(f"--{k}-xlsx"); ap.add_argument(f"--{k}-fasta")
    ap.add_argument("--curated-xlsx", help="Final_curated_SNPs_ONT_Illumina_v6.xlsx; adds the curated conditions")
    ap.add_argument("--curated-only", action="store_true", help="only the curated conditions (no raw tables, FASTA or BAMs needed)")
    ap.add_argument("--samples", help="samples.csv; curated isolate names are checked against it")
    ap.add_argument("--ill-bams", default=f"{B}/illumina_15x/step1/*/alignment_NC_002945v4/*_nodup.bam")
    ap.add_argument("--ont-bams", default=f"{B}/ont_15x/step1/*/alignment_NC_002945v4/*_nodup.bam")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chr", default="NC_002945.4")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-depth", type=float, default=5)
    ap.add_argument("--max-lowmq", type=float, default=0.25)
    ap.add_argument("--close", type=int, default=10, help="'close pair' = Illumina distance <= this")
    ap.add_argument("--expect", type=int, default=21)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.curated_only:
        if not a.curated_xlsx:
            sys.exit("--curated-only needs --curated-xlsx")
    elif not all(getattr(a, x) for x in ("ill_xlsx", "ill_fasta", "ont_xlsx", "ont_fasta")):
        sys.exit("need --ill-xlsx --ill-fasta --ont-xlsx --ont-fasta (or use --curated-only)")

    summ, tabs = [], {}
    if not a.curated_only:
        summ, tabs, _ = run_raw(a)
    if a.curated_xlsx:
        for label, (pi, po) in curated_conditions(a).items():
            s, c = summarise(pi, po, label, a); summ.append(s); tabs[label] = c
    res = pd.DataFrame(summ)
    res.to_csv(f"{a.out}/comparison_summary.tsv", sep="\t", index=False)
    txt = res.T.to_string(header=False)
    print("\n" + txt)
    open(f"{a.out}/distance_summary.txt", "w").write(
        f"Pairs per condition: {res.pairs.iloc[0]} (expected {a.expect * (a.expect - 1) // 2})\n"
        f"Close pairs = Illumina distance <= {a.close}\n\n" + txt + "\n")
    titles = {"unfiltered": "(a) vSNP3 alignment, no extra masking", "PE_PPE_only": "(b) PE/PPE sites removed",
              "PE_PPE_plus_callability": "(c) PE/PPE + callability filter",
              "curated": "(d) curated matrices (own sites)", "curated_shared_sites": "(e) curated, shared sites",
              "curated_complete_case": "(f) curated, complete-case sites"}
    figure3(tabs, res, titles, a)
    print("\nwrote Figure3_pairwise_SNP_distance.(png|pdf) and tables to", a.out)


if __name__ == "__main__":
    main()
