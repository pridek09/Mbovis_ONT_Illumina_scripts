#!/usr/bin/env python3
"""
STEP 2 - build samples.csv (isolate,platform,bam,zc_vcf) from the vSNP3 step 1 output folders.

Stops with an error unless exactly --expect (default 21) isolates are present on BOTH platforms,
each folder holds exactly one *_zc.vcf and one *_nodup.bam (no silent fallback to other BAMs).

Usage
  python 02_make_samples.py \
     --illumina /home/pridek/Illumina_ONT/analyses/vsnp3_supp/illumina_15x/step1 \
     --ont      /home/pridek/Illumina_ONT/analyses/vsnp3_supp/ont_15x/step1 \
     --out      samples.csv
"""
import argparse, csv, glob, os, re, sys

B = "/home/pridek/Illumina_ONT/analyses/vsnp3_supp"


def isolate_name(sample):
    s = re.sub(r"^\d{2}-\d{6}-", "", sample)                       # sequencing-run prefix
    s = re.sub(r"(_nanopore|_ont|_15x)+$", "", s, flags=re.I)      # platform / depth suffix
    m = re.match(r"^(ERR\d+)(?:ERR\d+)+$", s)                      # ONT runs with two glued accessions -> first one
    return m.group(1) if m else s


def scan(platform_dir):
    if not os.path.isdir(platform_dir):
        sys.exit(f"ERROR: {platform_dir} not found")
    found, problems = {}, []
    for d in sorted(os.listdir(platform_dir)):
        folder = os.path.join(platform_dir, d)
        if not os.path.isdir(folder) or d in ("Position_catalogs", "SNP_catalogs"):
            continue
        vcfs = sorted(glob.glob(os.path.join(folder, "**", "*_zc.vcf"), recursive=True))
        bams = sorted(glob.glob(os.path.join(folder, "**", "*_nodup.bam"), recursive=True))
        if len(vcfs) != 1 or len(bams) != 1:
            problems.append(f"{d}: {len(vcfs)} *_zc.vcf, {len(bams)} *_nodup.bam (need exactly 1 of each)")
            continue
        iso = isolate_name(d)
        if iso in found:
            problems.append(f"duplicate isolate name {iso} ({d} and {found[iso][0]})")
            continue
        found[iso] = (d, os.path.abspath(bams[0]), os.path.abspath(vcfs[0]))
    return found, problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--illumina", default=f"{B}/illumina_15x/step1")
    ap.add_argument("--ont", default=f"{B}/ont_15x/step1")
    ap.add_argument("--out", default="samples.csv")
    ap.add_argument("--expect", type=int, default=21)
    a = ap.parse_args()

    ill, p1 = scan(a.illumina)
    ont, p2 = scan(a.ont)
    paired = sorted(set(ill) & set(ont))
    print(f"Illumina folders: {len(ill)} | ONT folders: {len(ont)} | paired isolates: {len(paired)}")
    for x in sorted(set(ill) - set(ont)): print("  only Illumina:", x)
    for x in sorted(set(ont) - set(ill)): print("  only ONT:", x)
    for p in p1 + p2: print("  PROBLEM:", p)
    if p1 or p2 or len(paired) != a.expect or len(ill) != a.expect or len(ont) != a.expect:
        sys.exit(f"ERROR: expected exactly {a.expect} isolates on both platforms with no problems; not writing {a.out}")

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["isolate", "platform", "bam", "zc_vcf"])
        for iso in paired:
            w.writerow([iso, "illumina", ill[iso][1], ill[iso][2]])
            w.writerow([iso, "ont", ont[iso][1], ont[iso][2]])
    print(f"PASS: {a.expect} paired isolates; wrote {a.out}")


if __name__ == "__main__":
    main()
