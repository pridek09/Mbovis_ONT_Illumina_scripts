#!/usr/bin/env python3
"""
STEP 3 - per-isolate metrics (Table 1) and Illumina-vs-ONT SNP concordance (Tables 2-4),
computed exactly as defined in the manuscript Methods. Standard library only.

Usage
  python 03_snp_metrics_and_concordance.py --samples samples.csv --fai NC_002945v4.fa.fai --out results/

samples.csv (header required), one row per isolate x platform:
  isolate,platform,bam,zc_vcf
'bam' may also be a text file already holding `samtools depth -a` output (.txt/.depth/.gz).
platform must be 'illumina' or 'ont'.

Definitions (from the Methods)
  high-confidence SNP : _zc.vcf record, single-base REF and ALT in ACGT, INFO token AC=2, QUAL>=20
  ambiguous record    : INFO contains the exact token AC=1 (any record type; AC=1,1 not counted)
  uncalled position   : zero-coverage record (REF=N, ALT='.')
  quality SNP (Tab.1) : single-base REF/ALT, QUAL>150, any AC
  class of the OTHER platform at a position, in this order:
        shared > uncalled > ambiguous > reference
        ("reference" = none of the above: no high-confidence SNP, no zero-coverage record and no AC=1 record;
         this includes positions with a lower-quality or non-SNP record)
  QUAL offset         : vSNP3 raises ONT QUAL by 100. --ont-qual-offset N subtracts N from ONT QUAL before
                        the QUAL thresholds are applied (default 0 = use QUAL as written in the VCF).
                        qual_diagnostics.csv lists the QUAL distribution per sample so the offset can be checked.
  SDs are sample SDs (n - 1). Stops unless exactly --expect (21) isolates have both platforms.
"""
import argparse, csv, gzip, math, os, statistics as st, subprocess, sys
from collections import defaultdict

BASES = set("ACGT")


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_vcf(path, qual_offset=0.0):
    hc, amb, zero = {}, set(), set()
    n_quality, n_amb_records = 0, 0
    quals = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, ref, alt, qual, info = f[0], int(f[1]), f[3], f[4], f[5], f[7]
            if ref == "N" and alt == ".":
                zero.add((chrom, pos))
                continue
            tokens = info.split(";")
            try:
                q = float(qual)
                q = q - qual_offset
            except (ValueError, TypeError):
                q = None
            if "AC=1" in tokens:
                amb.add((chrom, pos))
                n_amb_records += 1
            single = len(ref) == 1 and len(alt) == 1 and ref in BASES and alt in BASES
            if single and "AC=2" in tokens and q is not None:
                quals.append(q + qual_offset)   # raw QUAL as written in the VCF
            if single and q is not None and q > 150:
                n_quality += 1
            if single and "AC=2" in tokens and q is not None and q >= 20:
                hc[(chrom, pos)] = alt
    return dict(hc=hc, amb=amb, zero=zero, n_quality=n_quality, n_amb=n_amb_records, quals=quals)


def depth_stats(path, L):
    if path.endswith(".bam"):
        proc = subprocess.Popen(["samtools", "depth", "-a", path], stdout=subprocess.PIPE, text=True)
        stream = proc.stdout
    else:
        proc, stream = None, open_text(path)
    total = zero = n = 0
    for line in stream:
        d = int(line.rsplit("\t", 1)[1])
        total += d
        n += 1
        if d == 0:
            zero += 1
    if proc:
        proc.wait()
        if proc.returncode:
            sys.exit(f"samtools failed on {path}")
    if n != L:
        print(f"WARNING {path}: depth file has {n} positions, reference length is {L}", file=sys.stderr)
    return total / L, zero


def ref_length(fai):
    with open(fai) as fh:
        return sum(int(l.split("\t")[1]) for l in fh if l.strip())


def classify(pos, other):
    if pos in other["hc"]:
        return "shared"
    if pos in other["zero"]:
        return "uncalled"
    if pos in other["amb"]:
        return "ambiguous"
    return "reference"


def pct(a, b):
    return 100.0 * a / b if b else float("nan")


def msd(values):
    v = [x for x in values if not math.isnan(x)]
    return f"{st.mean(v):.1f} ± {st.stdev(v):.1f}" if len(v) > 1 else ""


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--fai", help="reference .fai (gives L); otherwise --length")
    ap.add_argument("--length", type=int, default=4349904, help="reference length if no .fai")
    ap.add_argument("--out", default="results")
    ap.add_argument("--expect", type=int, default=21)
    ap.add_argument("--ont-qual-offset", type=float, default=0.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    L = ref_length(a.fai) if a.fai else a.length
    print("reference length L =", L)

    data = defaultdict(dict)
    for r in csv.DictReader(open(a.samples)):
        plat = r["platform"].strip().lower()
        v = parse_vcf(r["zc_vcf"], a.ont_qual_offset if plat == "ont" else 0.0)
        v["avg_depth"], v["zero_cov"] = depth_stats(r["bam"], L)
        v["cov_pct"] = 100.0 * (L - v["zero_cov"]) / L
        v["zero_pct"] = 100.0 * v["zero_cov"] / L
        if len(v["zero"]) != v["zero_cov"]:
            print(f"NOTE {r['isolate']} {plat}: {len(v['zero'])} zc records vs "
                  f"{v['zero_cov']} zero-depth positions (check the _zc.vcf)", file=sys.stderr)
        data[r["isolate"]][plat] = v

    isolates = sorted(data)
    missing = [i for i in isolates if set(data[i]) != {"illumina", "ont"}]
    if missing:
        sys.exit(f"isolates without both platforms: {missing}")
    if len(isolates) != a.expect:
        sys.exit(f"ERROR: expected {a.expect} isolates, found {len(isolates)}")

    qrows = []
    for i in isolates:
        for p in ("illumina", "ont"):
            q = sorted(data[i][p]["quals"])
            if q:
                qrows.append([i, p, len(q), f"{q[0]:.1f}", f"{st.median(q):.1f}", f"{q[-1]:.1f}"])
    write_csv(f"{a.out}/qual_diagnostics.csv",
              ["isolate", "platform", "n_single_base_AC2_records", "QUAL_min", "QUAL_median", "QUAL_max"], qrows)
    om = [float(r[3]) for r in qrows if r[1] == "ont"]
    im = [float(r[3]) for r in qrows if r[1] == "illumina"]
    if om and im:
        print(f"QUAL minimum over AC=2 SNP records: Illumina {min(im):.0f}, ONT {min(om):.0f} "
              f"(offset applied to ONT: {a.ont_qual_offset:g}). If ONT minimum is ~100 higher than Illumina, "
              "rerun with --ont-qual-offset 100.")

    rows = []
    for i in isolates:
        for p in ("illumina", "ont"):
            v = data[i][p]
            rows.append([i, p, f"{v['avg_depth']:.2f}", f"{v['cov_pct']:.2f}", v["zero_cov"],
                         f"{v['zero_pct']:.3f}", len(v["hc"]), v["n_amb"], v["n_quality"]])
    write_csv(f"{a.out}/table1_per_isolate_metrics.csv",
              ["isolate", "platform", "avg_depth", "genome_with_coverage_pct", "no_coverage_bases",
               "no_coverage_pct", "high_conf_snps", "ambiguous_snps_AC1", "quality_snps_QUAL_gt150"], rows)

    t2, t3, t4, pos_rows = [], [], [], []
    tot = defaultdict(int)
    pc = defaultdict(list)
    for i in isolates:
        ill, ont = data[i]["illumina"], data[i]["ont"]
        sh = set(ill["hc"]) & set(ont["hc"])
        ill_only = set(ill["hc"]) - set(ont["hc"])
        ont_only = set(ont["hc"]) - set(ill["hc"])
        same = sum(1 for p in sh if ill["hc"][p] == ont["hc"][p])

        ic = defaultdict(int)
        oc = defaultdict(int)
        for p in sorted(ill_only):
            c = classify(p, ont); ic[c] += 1
            pos_rows.append([i, p[0], p[1], "illumina_only", c, ill["hc"][p], ""])
        for p in sorted(ont_only):
            c = classify(p, ill); oc[c] += 1
            pos_rows.append([i, p[0], p[1], "ont_only", c, "", ont["hc"][p]])
        for p in sorted(sh):
            pos_rows.append([i, p[0], p[1], "shared", "shared", ill["hc"][p], ont["hc"][p]])

        ill_val, ont_val = ic["reference"], oc["reference"]
        excluded = (ic["ambiguous"] + ic["uncalled"]) + (oc["ambiguous"] + oc["uncalled"])
        validated = len(sh) + ill_val + ont_val
        allpos = len(sh) + len(ill_only) + len(ont_only)
        conc, conc_noex = pct(len(sh), validated), pct(len(sh), allpos)
        rec, prec = pct(len(sh), len(ill["hc"])), pct(len(sh), len(ont["hc"]))

        t2.append([i, len(ill["hc"]), len(ont["hc"]), len(sh), same, f"{pct(same, len(sh)):.1f}"])
        t3.append([i, len(sh), ill_val, ont_val, excluded, validated, f"{conc:.1f}", f"{conc_noex:.1f}"])
        t4.append([i, len(ill["hc"]), len(ont["hc"]), len(sh), len(ill_only), len(ont_only),
                   f"{rec:.1f}", f"{prec:.1f}", ic["reference"], ic["ambiguous"], ic["uncalled"],
                   oc["reference"], oc["ambiguous"], oc["uncalled"]])
        for k, val in dict(ill=len(ill["hc"]), ont=len(ont["hc"]), shared=len(sh), same=same,
                           ill_val=ill_val, ont_val=ont_val, excluded=excluded, validated=validated,
                           allpos=allpos, ill_only=len(ill_only), ont_only=len(ont_only),
                           ic_ref=ic["reference"], ic_amb=ic["ambiguous"], ic_unc=ic["uncalled"],
                           oc_ref=oc["reference"], oc_amb=oc["ambiguous"], oc_unc=oc["uncalled"]).items():
            tot[k] += val
        pc["allelic"].append(pct(same, len(sh)))
        pc["conc"].append(conc); pc["conc_noex"].append(conc_noex)
        pc["rec"].append(rec); pc["prec"].append(prec)

    T = tot
    t2.append(["Pooled", T["ill"], T["ont"], T["shared"], T["same"], f"{pct(T['same'], T['shared']):.1f}"])
    t2.append(["Mean ± SD", "", "", "", "", msd(pc["allelic"])])
    t3.append(["Pooled", T["shared"], T["ill_val"], T["ont_val"], T["excluded"], T["validated"],
               f"{pct(T['shared'], T['validated']):.1f}", f"{pct(T['shared'], T['allpos']):.1f}"])
    t3.append(["Mean ± SD", "", "", "", "", "", msd(pc["conc"]), msd(pc["conc_noex"])])
    t4.append(["Pooled", T["ill"], T["ont"], T["shared"], T["ill_only"], T["ont_only"],
               f"{pct(T['shared'], T['ill']):.1f}", f"{pct(T['shared'], T['ont']):.1f}",
               T["ic_ref"], T["ic_amb"], T["ic_unc"], T["oc_ref"], T["oc_amb"], T["oc_unc"]])
    t4.append(["Mean ± SD", "", "", "", "", "", msd(pc["rec"]), msd(pc["prec"])] + [""] * 6)

    write_csv(f"{a.out}/table2_allelic_concordance.csv",
              ["isolate", "illumina_hc_snps", "ont_hc_snps", "shared", "same_alt", "allelic_concordance_pct"], t2)
    write_csv(f"{a.out}/table3_position_concordance.csv",
              ["isolate", "shared", "illumina_only_validated", "ont_only_validated", "excluded",
               "validated", "position_concordance_pct", "concordance_without_exclusions_pct"], t3)
    write_csv(f"{a.out}/table4_recovery_precision.csv",
              ["isolate", "illumina_snps", "ont_snps", "shared", "illumina_only", "ont_only",
               "recovery_pct", "ont_precision_pct", "ill_only_ONT_ref", "ill_only_ONT_ambig",
               "ill_only_ONT_uncalled", "ont_only_ILL_ref", "ont_only_ILL_ambig", "ont_only_ILL_uncalled"], t4)
    write_csv(f"{a.out}/snp_positions_all_classes.csv",
              ["isolate", "contig", "pos", "set", "other_platform_call", "illumina_alt", "ont_alt"], pos_rows)

    print("\nPooled check:  Illumina", T["ill"], "| ONT", T["ont"], "| shared", T["shared"],
          "| validated", T["validated"], "| excluded", T["excluded"],
          f"| position concordance {pct(T['shared'], T['validated']):.1f}% "
          f"| without exclusions {pct(T['shared'], T['allpos']):.1f}%")


if __name__ == "__main__":
    main()
