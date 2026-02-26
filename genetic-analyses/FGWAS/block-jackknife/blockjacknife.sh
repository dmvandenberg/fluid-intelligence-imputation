#!/usr/bin/env bash
#SBATCH --job-name=jk_corr
#SBATCH --array=1-200
#SBATCH --partition=rome
#SBATCH --time=00:30:00
#SBATCH --output=/home/dvandenberg/project_folder/iq-from-scratch/downstream/F-GWAS/jackknife/logs/jk_%a.out

set -euo pipefail

# ---- CONFIG ----
SNIPAR_PY="/home/dvandenberg/storage/software/snipar-env-patch/bin/python"
CORRELATE="/home/dvandenberg/storage/software/snipar-env-patch/bin/correlate.py"

SUMSTATS="/home/dvandenberg/project_folder/iq-from-scratch/downstream/F-GWAS/raw_sumstats/average_mega.sumstats.gz"
QC_SNPLIST="/home/dvandenberg/project_folder/iq-from-scratch/downstream/F-GWAS/snps.098.ids"
BLOCK_DIR="/home/dvandenberg/project_folder/iq-from-scratch/downstream/F-GWAS/jackknife/blocks"

# correlate.py arguments
OUT_BASE="/home/dvandenberg/project_folder/iq-from-scratch/downstream/F-GWAS/jackknife/average_mega"
LD_PREFIX="/home/dvandenberg/project_folder/iq-from-scratch/downstream/F-GWAS/ldscores/chr_@"
THREADS=16

# ----- RUN -----
N="${SLURM_ARRAY_TASK_ID}"
EXCL="${BLOCK_DIR}/block_${N}.snps"

WORK="${TMPDIR:-/tmp}/jk_${N}"
mkdir -p "${WORK}"

FILTERED="${WORK}/minus_block_${N}.sumstats.gz"
FILTERED_PREFIX="${WORK}/minus_block_${N}"

echo "block ${N}"
echo "QC list: ${QC_SNPLIST}"
echo "Excluding SNPs listed in: ${EXCL}"

# Keep header; exclude rows where SNP (col2) is in EXCL
zcat "${SUMSTATS}" | awk -v excl="${EXCL}" -v qclist="${QC_SNPLIST}" '
  BEGIN {
    while ((getline s < qclist) > 0) qc[s]=1
    close(qclist)
    while ((getline s < excl) > 0) ex[s]=1
    close(excl)
  }
  NR==1 { print; next }
  (qc[$2] && !ex[$2]) { print }
' | gzip -c > "${FILTERED}"

KEPT=$(zcat "${FILTERED}" | awk "NR>1{c++} END{print c+0}")
echo "SNPs kept after QC+exclude: ${KEPT}"

# Output per block
OUT_PREFIX="${OUT_BASE}/block_${N}"

export NUMBA_DISABLE_JIT=1
"${SNIPAR_PY}" "${CORRELATE}" \
  "${FILTERED_PREFIX}" \
  "${OUT_PREFIX}" \
  --ldscores "${LD_PREFIX}" \
  --threads "${THREADS}"

echo "Done block ${N}. Output: ${OUT_PREFIX}"

rm -rf "${WORK}"