#!/usr/bin/env bash
#SBATCH --job-name=mk_blocks
#SBATCH --partition=rome
#SBATCH --time=00:30:00
#SBATCH --output=../F-GWAS/jackknife/logs/mk_blocks.out

#!/usr/bin/env bash
set -euo pipefail

# ---- CONFIG ----
GWAS1="../F-GWAS/raw_sumstats/average_mega.sumstats.gz"
GWAS2="../F-GWAS/raw_sumstats/average_measured.sumstats.gz"
GWAS3="../F-GWAS/raw_sumstats/fis1_measured.sumstats.gz"

OUT_BLOCK_DIR="../F-GWAS/jackknife/blocks"
QC_SNPLIST="../F-GWAS/snps.098.ids"
N_BLOCKS=200

# ---- SETUP ----
TMPDIR="$(mktemp -d)"
trap 'rm -rf "${TMPDIR}"' EXIT

OUT_SNPS="${OUT_BLOCK_DIR}/all_snps.txt"
INTERSECT_SNPS="${TMPDIR}/intersect.snps"

echo "Building SNP intersection across 3 GWAS, constrained to QC list..."

# Helper: extract SNP IDs (col2) present in sumstats AND in QC list; output sorted unique SNPs
extract_snps_sorted_unique () {
  local in_gwas="$1"
  local out_txt="$2"

  zcat "${in_gwas}" \
    | awk -v qclist="${QC_SNPLIST}" '
        BEGIN {
          while ((getline s < qclist) > 0) qc[s]=1
          close(qclist)
        }
        NR==1 { next }          # skip header
        qc[$2] { print $2 }     # SNP id in col2
      ' \
    | sort -u > "${out_txt}"
}

extract_snps_sorted_unique "${GWAS1}" "${TMPDIR}/gwas1.snps"
extract_snps_sorted_unique "${GWAS2}" "${TMPDIR}/gwas2.snps"
extract_snps_sorted_unique "${GWAS3}" "${TMPDIR}/gwas3.snps"

# Intersection on SNP IDs
comm -12 "${TMPDIR}/gwas1.snps" "${TMPDIR}/gwas2.snps" \
  | comm -12 - "${TMPDIR}/gwas3.snps" \
  > "${INTERSECT_SNPS}"

N_INTERSECT=$(wc -l < "${INTERSECT_SNPS}")
echo "Intersect SNPs (QC ∩ all 3 GWAS): ${N_INTERSECT}"

if [[ "${N_INTERSECT}" -eq 0 ]]; then
  echo "ERROR: Intersection is empty. Check SNP column (expected col2) and/or QC list compatibility." >&2
  exit 1
fi

echo "Ordering intersect SNPs by chr/pos using GWAS1 as anchor (${GWAS1})..."

# Create ordered SNP list by chr/pos from anchor GWAS1, filtered to intersect set
zcat "${GWAS1}" \
  | awk -v keep="${INTERSECT_SNPS}" '
      BEGIN {
        while ((getline s < keep) > 0) ok[s]=1
        close(keep)
      }
      NR==1 { next }
      ok[$2] { print $1, $3, $2 }   # CHR POS SNP
    ' \
  | sort -k1,1n -k2,2n \
  | awk '{ print $3 }' \
  > "${OUT_SNPS}"

N_SNPS=$(wc -l < "${OUT_SNPS}")
echo "Reference SNPs written (ordered): ${N_SNPS}"
echo "Output: ${OUT_SNPS}"

# ---- BLOCKS ----
BLOCK_SIZE=$(( (N_SNPS + N_BLOCKS - 1) / N_BLOCKS ))
echo "Splitting into ${N_BLOCKS} blocks of up to ${BLOCK_SIZE} SNPs..."

# Clean any old blocks matching pattern to avoid mixing (optional, but often desirable)
rm -f "${OUT_BLOCK_DIR}/block_"*.snps "${OUT_BLOCK_DIR}/block_"* 2>/dev/null || true

split -l "${BLOCK_SIZE}" -d -a 3 "${OUT_SNPS}" "${OUT_BLOCK_DIR}/block_"

# Rename block_000 -> block_1.snps ... sequentially
i=1
for f in "${OUT_BLOCK_DIR}"/block_*; do
  mv "${f}" "${OUT_BLOCK_DIR}/block_${i}.snps"
  i=$((i+1))
done

echo "Blocks created: $((i-1))"