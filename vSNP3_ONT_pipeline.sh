#!/bin/bash
# vSNP3 on the 15x ONT reads (single-end), following Supplementary Material S1
# and Ghielmetti et al. 2023 (Front. Microbiol. 14:1307440):
#   Step 1: minimap2 mapping, SAMtools fixmate/sort/markdup, BCFtools SNP calling (--nanopore);
#           vSNP3 raises ONT QUAL values by 100 so Illumina and ONT VCFs can be run together
#   Step 2: vSNP3 step 2 with -a (default filters; PE/PPE sites are masked afterwards in 05_alignment_distances.py), SNP table, SNP alignment, RAxML tree
#   Then:   pairwise SNP distances with snp-dists
# All vSNP3 options other than the documented ones are left at their defaults.
set -uo pipefail

# ------------------------------------------------------------------ settings
source "$(conda info --base)/etc/profile.d/conda.sh"
VSNP_ENV="vsnp3_clean"
SNPDIST_ENV="snp_dists_env"          # snp-dists (the paper used v0.8.2)
REFTYPE="Mycobacterium_AF2122"
RAW="/home/pridek/Illumina_ONT/downsampled/ONT/15x/ONT_15x"             # <sample>.fastq.gz (e.g. RAD24002_15x.fastq.gz)
OUT="/home/pridek/Illumina_ONT/analyses/vsnp3_supp/ont_15x"             # new folder: existing results are not touched
MIN_MB=10                                                               # skip tiny files
STEP1="$OUT/step1"
STEP2="$OUT/step2"
VCF_DIR="$STEP2/vcfs_flat"
# ---------------------------------------------------------------------------

conda activate "$VSNP_ENV" || { echo "cannot activate $VSNP_ENV"; exit 1; }
mkdir -p "$STEP1" "$STEP2"

echo "Checking reference type $REFTYPE"
if ! vsnp3_path_adder.py -s | grep -q "$REFTYPE"; then
    echo "ERROR: reference type '$REFTYPE' not found. Install it with:"
    echo "  vsnp3_path_adder.py -d /home/pridek/vsnp3_test_dataset/vsnp_dependencies"
    exit 1
fi

{
  echo "run date: $(date -Iseconds)"; echo "conda env: $CONDA_DEFAULT_ENV"; echo "reference type: $REFTYPE"
  for t in bwa samtools freebayes minimap2 bcftools raxmlHPC; do
      printf "%s: " "$t"; ( $t --version 2>&1 || $t version 2>&1 ) | head -1
  done
} > "$OUT/run_parameters.txt" 2>&1

# ---------------------------------------------------------------- Step 1
FAILED=()
for FASTQ in "$RAW"/*.fastq.gz; do
    SAMPLE=$(basename "$FASTQ" .fastq.gz)
    SIZE_MB=$(du -m "$FASTQ" | cut -f1)
    if [ "$SIZE_MB" -lt "$MIN_MB" ]; then
        echo "Skipping $SAMPLE (only ${SIZE_MB} MB)"; continue
    fi
    SAMPLE_OUT="$STEP1/$SAMPLE"
    if ls "$SAMPLE_OUT"/alignment_*/*_zc.vcf >/dev/null 2>&1; then
        echo "Skipping $SAMPLE (step 1 already done)"; continue
    fi
    mkdir -p "$SAMPLE_OUT"
    echo "===================================="
    echo "Step 1: $SAMPLE ($FASTQ)"
    echo "===================================="
    vsnp3_step1.py -n "$SAMPLE" -r1 "$FASTQ" -t "$REFTYPE" --nanopore -o "$SAMPLE_OUT" \
        > "$SAMPLE_OUT/vsnp3_step1.log" 2>&1 || { echo "FAILED: $SAMPLE"; FAILED+=("$SAMPLE"); }
done
if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "Step 1 failures: ${FAILED[*]}"
    echo "See vsnp3_step1.log in each sample folder. Fix them and re-run (finished samples are skipped)."
    exit 1
fi

# ---------------------------------------------------------------- Step 2
echo "===== vSNP3 STEP 2 ====="
rm -rf "$VCF_DIR"; mkdir -p "$VCF_DIR"
find "$STEP1" -path "*/alignment_*/*_zc.vcf" -exec cp {} "$VCF_DIR" \;
VCF_COUNT=$(find "$VCF_DIR" -name "*.vcf" | wc -l)
echo "Flattened VCF count: $VCF_COUNT"
[ "$VCF_COUNT" -ge 2 ] || { echo "ERROR: fewer than 2 VCFs"; exit 1; }

vsnp3_step2.py -a -t "$REFTYPE" -wd "$VCF_DIR" -o "$STEP2" 2>&1 | tee "$STEP2/vsnp3_step2.log"

# ---------------------------------------------------------------- SNP distances
conda deactivate
conda activate "$SNPDIST_ENV" || { echo "cannot activate $SNPDIST_ENV"; exit 1; }
echo "snp-dists version: $(snp-dists -v 2>&1)" >> "$OUT/run_parameters.txt"
ALN=$(find "$STEP2" -name "*.fasta" -printf '%T@ %p\n' | grep -i "Mbovis-All\|-All-" | sort -n | tail -1 | cut -d' ' -f2-)
[ -n "$ALN" ] || { echo "ERROR: SNP alignment fasta not found under $STEP2"; exit 1; }
echo "SNP alignment: $ALN"
snp-dists "$ALN" > "${ALN%.fasta}_snp_distances.tsv"
echo "Pairwise SNP distances: ${ALN%.fasta}_snp_distances.tsv"
echo "===== DONE ====="
