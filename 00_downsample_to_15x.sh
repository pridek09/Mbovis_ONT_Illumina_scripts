#!/usr/bin/env bash
# STEP 0 - downsample all Illumina and ONT datasets to 15x with Rasusa (seed 42).
# Mbovis AF2122/97 (NC_002945.4) genome length = 4,349,904 bp.
#
# Writes the downsampled reads to the folders that vSNP3_pipeline.sh and vSNP3_ONT.sh read from,
# so the reads analysed by vSNP3 are exactly the reads this script produced.
# Everything needed for the Methods is logged: Rasusa version, seed, per-sample Rasusa logs,
# and a seqkit table with the realised depth (total bases / genome size) of input and output.
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"

BASE=/home/pridek/Illumina_ONT
RAW=$BASE/raw_data
COV=15
GENOME=4349904
SEED=42

ILL_IN="$RAW/Illumina/illumina"          # <sample>_1.fastq.gz / <sample>_2.fastq.gz
ONT_IN="$RAW/ONT"                        # <sample>.fastq.gz

ILL_OUT="$BASE/downsampled/Illumina/15x/illumina_15x"   # read by vSNP3_pipeline.sh
ONT_OUT="$BASE/downsampled/ONT/15x/ONT_15x"             # read by vSNP3_ONT.sh
LOGDIR="$BASE/downsampled/logs"
mkdir -p "$ILL_OUT" "$ONT_OUT" "$LOGDIR"

conda activate rasusa_env
{
  echo "date: $(date -Iseconds)"
  echo "rasusa: $(rasusa --version)"
  echo "coverage: $COV"
  echo "genome_size: $GENOME"
  echo "seed: $SEED"
} | tee "$LOGDIR/rasusa_run_parameters.txt"

# ---------------- Illumina (paired-end)
for R1 in "$ILL_IN"/*_1.fastq.gz; do
    S=$(basename "$R1" _1.fastq.gz)
    echo "Illumina $S"
    rasusa reads --coverage "$COV" --genome-size "$GENOME" --seed "$SEED" \
        "$R1" "$ILL_IN/${S}_2.fastq.gz" \
        -o "$ILL_OUT/${S}_1_15x.fastq.gz" -o "$ILL_OUT/${S}_2_15x.fastq.gz" \
        2> "$LOGDIR/${S}.illumina.rasusa.log"
done

# ---------------- ONT (single-end)
for FQ in "$ONT_IN"/*.fastq.gz; do
    S=$(basename "$FQ" .fastq.gz)
    echo "ONT $S"
    rasusa reads --coverage "$COV" --genome-size "$GENOME" --seed "$SEED" \
        "$FQ" -o "$ONT_OUT/${S}_15x.fastq.gz" \
        2> "$LOGDIR/${S}.ont.rasusa.log"
done

# ---------------- statistics of inputs and outputs
conda activate seqkit-env
seqkit stats -T "$ILL_IN"/*.fastq.gz  > "$LOGDIR/illumina_input_seqkit_stats.tsv"
seqkit stats -T "$ONT_IN"/*.fastq.gz  > "$LOGDIR/ont_input_seqkit_stats.tsv"
seqkit stats -T "$ILL_OUT"/*.fastq.gz > "$LOGDIR/illumina_15x_seqkit_stats.tsv"
seqkit stats -T "$ONT_OUT"/*.fastq.gz > "$LOGDIR/ont_15x_seqkit_stats.tsv"

# realised depth = sum_bases / genome size (Illumina: R1 + R2 of a sample summed)
python3 - "$GENOME" "$LOGDIR" <<'PY'
import csv, re, sys, collections
G = int(sys.argv[1]); d = sys.argv[2]
rows = []
for tag in ("illumina_input", "ont_input", "illumina_15x", "ont_15x"):
    tot = collections.defaultdict(int)
    for r in csv.DictReader(open(f"{d}/{tag}_seqkit_stats.tsv"), delimiter="\t"):
        s = re.sub(r"_[12](_15x)?\.fastq\.gz$|(_15x)?\.fastq\.gz$", "", r["file"].split("/")[-1])
        tot[s] += int(r["sum_len"].replace(",", ""))
    rows += [(tag, s, b, round(b / G, 2)) for s, b in sorted(tot.items())]
with open(f"{d}/realised_depth.tsv", "w") as fh:
    fh.write("set\tsample\ttotal_bases\tdepth_x\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
low = [r for r in rows if r[0].endswith("_input") and r[3] < 15]
print(f"samples with input depth < 15x (kept in full by Rasusa): {len(low)}")
for r in low: print("  ", r)
print("wrote", f"{d}/realised_depth.tsv")
PY
echo "Finished. Reads: $ILL_OUT  $ONT_OUT   Logs: $LOGDIR"
