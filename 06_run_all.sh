#!/usr/bin/env bash
# Run the whole analysis in order. Every step stops on error. Re-running overwrites results.
#   bash 06_run_all.sh            # steps 2-5 (analysis only; vSNP3 runs are expensive)
#   bash 06_run_all.sh --full     # also runs downsampling (00) and both vSNP3 pipelines (01)
# Curated matrices: put Final_curated_SNPs_ONT_Illumina_v6.xlsx in $BASE, or set CURATED_XLSX=/path/to/file.xlsx
set -euo pipefail
S="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE=/home/pridek/Illumina_ONT
A=$BASE/analyses/vsnp3_supp
RES=$A/results_final
mkdir -p "$RES"

if [[ "${1:-}" == "--full" ]]; then
  bash "$S/00_downsample_to_15x.sh"
  bash "$S/vSNP3_pipeline.sh"
  bash "$S/vSNP3_ONT.sh"
fi
bash "$S/01_record_versions.sh"

python3 "$S/02_make_samples.py" --out "$RES/samples.csv"
python3 "$S/03_snp_metrics_and_concordance.py" --samples "$RES/samples.csv" --ont-qual-offset 100 \
    --out "$RES/concordance" 2>&1 | tee "$RES/concordance.log"
python3 "$S/04_stats_and_figures.py" --results "$RES/concordance" --out "$RES/stats" | tee "$RES/stats.log"

CUR=${CURATED_XLSX:-$BASE/Final_curated_SNPs_ONT_Illumina_v6.xlsx}   # curated matrices (sheets Cleaned_ONT / Cleaned_Illumina)
CURARGS=()
if [[ -f "$CUR" ]]; then CURARGS=(--curated-xlsx "$CUR" --samples "$RES/samples.csv"); else echo "NOTE: $CUR not found - curated conditions skipped"; fi

latest() { ls -t "$1"/$2 2>/dev/null | head -n1; }
ILL=$A/illumina_15x/step2/Mbovis-All; ONT=$A/ont_15x/step2/Mbovis-All
python3 "$S/05_alignment_distances.py" \
  --ill-xlsx  "$(latest "$ILL" '*sorted_table*.xlsx')" --ill-fasta "$(latest "$ILL" 'Mbovis-All-*.fasta')" \
  --ont-xlsx  "$(latest "$ONT" '*sorted_table*.xlsx')" --ont-fasta "$(latest "$ONT" 'Mbovis-All-*.fasta')" \
  "${CURARGS[@]}" --out "$RES/distances" | tee "$RES/distances.log"
echo "All done. Results in $RES. Read concordance.log, stats.log, distances.log for warnings."
