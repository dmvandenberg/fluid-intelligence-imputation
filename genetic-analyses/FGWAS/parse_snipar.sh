#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --partition=rome
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --output=./mean.iq/parse.log

module load 2025
module load R/4.4.2-gfbf-2025a

# Locations of SNIPAR python version and correlate.py
SNIPAR_PY="./snipar-env-patch/bin/python"
CORRELATE="./snipar-env-patch/bin/correlate.py"
LD_PREFIX="./eu_ld/_chr_@"

# Change this to folder where all per chromosome sumstats are downloaded 
# from UKB RAP
TRAIT="mean.iq"
IN_DIR="./${TRAIT}"

# Merge chromosomes into one file
cd ${IN_DIR}
OUT="${TRAIT}.sumstats.gz"
zcat "${IN_DIR}/chr1.sumstats.gz" | head -n 1 | gzip > "${OUT}"

for chr in {1..22}; do
  zcat "${IN_DIR}/chr${chr}.sumstats.gz" | tail -n +2
done | gzip >> "${OUT}"

# Check N snps
zcat ${TRAIT}.sumstats.gz | wc -l

# Run correlate.py
export NUMBA_DISABLE_JIT=1
"${SNIPAR_PY}" "${CORRELATE}" \
  ${TRAIT} \
  "${IN_DIR}/${TRAIT}" \
  --ldscores "${LD_PREFIX}" \
  --threads "${THREADS}"

IN_F="$(realpath "${IN_DIR}/${OUT}")"

# Split into direct_effect and NTC sumstats
Rscript --vanilla - ${IN_F} ${IN_DIR} ${TRAIT} <<'EOF'
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
f <- args[1]
p <- args[2]
t <- args[3]

con <- gzfile(f, open = "rt")
file <- read.table(con, header = TRUE, stringsAsFactors = FALSE)
close(con)

NTC <- file %>%
  select(SNP, CHR=chromosome, POS=pos, A1, A2, FRQ=freq,
         BETA=avg_NTC_Beta, SE=avg_NTC_SE, Z=avg_NTC_Z, N=avg_NTC_N) %>%
  mutate(P = 2 * pnorm(-abs(Z))) %>%
  select(-Z)

direct <- file %>%
  select(SNP, CHR=chromosome, POS=pos, A1, A2, FRQ=freq,
         BETA=direct_Beta, SE=direct_SE, Z=direct_Z, N=direct_N) %>%
  mutate(P = 2 * pnorm(-abs(Z))) %>%
  select(-Z)

write.table(NTC, paste0(p,"/", t, ".NTC.gwas.txt"), sep="\t", row.names=FALSE, quote=FALSE)
write.table(direct, paste0(p,"/", t, ".direct.gwas.txt"), sep="\t", row.names=FALSE, quote=FALSE)
EOF
