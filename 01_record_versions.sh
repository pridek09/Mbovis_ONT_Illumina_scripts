#!/usr/bin/env bash
# STEP 1b - record every software version for the Methods (run AFTER the two vSNP3 pipelines).
# The old run_parameters.txt files captured error text for bwa and raxmlHPC and never recorded vSNP3.
# Output: analyses/vsnp3_supp/software_versions.txt
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
OUT=/home/pridek/Illumina_ONT/analyses/vsnp3_supp/software_versions.txt
mkdir -p "$(dirname "$OUT")"

first() { "$@" 2>&1 | head -n1; }
{
  echo "date: $(date -Iseconds)"

  conda activate vsnp3_clean
  echo "== conda env vsnp3_clean =="
  conda list -n vsnp3_clean 2>/dev/null | grep -Ei '^(vsnp3|bwa|samtools|bcftools|freebayes|minimap2|raxml|snp-dists|biopython|pysam) ' || true
  echo "bwa:        $(bwa 2>&1 | grep -i '^Version' | head -n1)"
  echo "samtools:   $(first samtools --version)"
  echo "bcftools:   $(first bcftools --version)"
  echo "freebayes:  $(first freebayes --version)"
  echo "minimap2:   $(first minimap2 --version)"
  echo "raxmlHPC:   $(raxmlHPC -v 2>&1 | grep -i 'version' | head -n1)"
  echo "vsnp3_step1: $(first vsnp3_step1.py --version)"
  conda deactivate

  conda activate snp_dists_env
  echo "snp-dists:  $(first snp-dists -v)"
  conda deactivate

  conda activate rasusa_env
  echo "rasusa:     $(first rasusa --version)"
  conda deactivate

  conda activate seqkit-env
  echo "seqkit:     $(first seqkit version)"
  conda deactivate

  echo "== analysis environment =="
  python3 - <<'PY'
import sys, scipy, numpy, pandas, openpyxl, matplotlib
print("python", sys.version.split()[0], "| scipy", scipy.__version__, "| numpy", numpy.__version__,
      "| pandas", pandas.__version__, "| openpyxl", openpyxl.__version__, "| matplotlib", matplotlib.__version__)
PY
  echo "samtools (analysis): $(first samtools --version)"
} | tee "$OUT"
echo "wrote $OUT  -- check that every line has a version, and copy them into the Methods"
