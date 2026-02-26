#!/usr/bin/env bash
set -euo pipefail

###############################################################
# Run Regenie GWAS STEP 1 + 2 (robust SEs, dummy GxE)
###############################################################

# (Optional) modules on Snellius
module load 2025 || true
module load Python/3.13.1-GCCcore-14.2.0 || true
module load SciPy-bundle/2025.06-gfbf-2025a || true

###############################################################
# User input
###############################################################

# Project (NOTE trailing colon)
PROJECT="project-J34ZPG8JK1BB2FFyy11gB5Qk:"

# Input files
COV="${PROJECT}/input/fis.cov"
PHENO="${PROJECT}/input/fis.all.complete.pheno"
SNP_QC_STEP1="${PROJECT}/genotypes/qc_pass.snplist"
SNP_QC_STEP2="${PROJECT}/genotypes/HRC/ukb_imp.HRC.v2.snplist"
GENO_PREFIX="${PROJECT}/genotypes/merged_noQC/ukb22418_c1_22_v2_merged"
IMP_PREFIX="${PROJECT}/Bulk/Imputation/UKB imputation from genotype"

# Output directories
STEP1_OUT="${PROJECT}/output/regenie_step_1_all_complete/"
STEP1_OUT_KEEP="${PROJECT}/output/step1_keep.txt"
STEP2_OUT="${PROJECT}/output/regenie_step_2_HRC_all_complete/"

# Run on
INSTANCE="mem1_ssd1_v2_x16"

# OPTIONAL IF STEP1 DONE WITH SAME PHENOTYPES + INDIVIDUALS 
STEP1=$(dx run swiss-army-knife \
  -imount_inputs=true \
  -iin="${GENO_PREFIX}.bed" \
  -iin="${GENO_PREFIX}.bim" \
  -iin="${GENO_PREFIX}.fam" \
  -iin="${PHENO}" \
  -iin="${COV}" \
  -iin="${SNP_QC_STEP1}" \
  -icmd="bash -lc 'set -euo pipefail; regenie --step 1 --bed ukb22418_c1_22_v2_merged --phenoFile fis.all.complete.pheno --extract qc_pass.snplist --covarFile fis.cov --phenoColList imputed.fis2,fis2,all.fis2,iq2,imputed.iq2,all.iq2,earliest.iq,imputed.earliest.iq,all.earliest.iq2,mean.iq,imputed.mean.iq,all.mean.iq2,mean.iq.scaled,MEGAV1,MEGAV2,impvA,impvB,all.iq2.never,all.iq2.ever --covarColList pca1,pca2,pca3,pca4,pca5,pca6,pca7,pca8,pca9,pca10,pca11,pca12,pca13,pca14,pca15,pca16,pca17,pca18,pca19,pca20,pca21,pca22,pca23,pca24,pca25,age2,age.squared2,sex2,agebysex2,age.squaredbysex2,array,E1 --bsize 1000 --threads 48 --lowmem --gz --out fis_step_1'" \
  --name "Regenie Step 1 HRC" \
  --instance-type "${INSTANCE}" \
  --destination "${STEP1_OUT}" \
  --priority high \
  --brief -y)

dx wait "${STEP1}"
echo "Step 1 complete."

# Run step 2. 

# Helper for generating download commpand for .loco files from step1
DOWNLOAD=""
for i in {1..19}; do
DOWNLOAD+=" dx download ${STEP1_OUT}fis_step_1_${i}.loco.gz -o /home/dnanexus/out/out/;"
done

for chr in {1..22}; do
STEP2=$(dx run swiss-army-knife \
  -imount_inputs=true \
  -iin="${IMP_PREFIX}/ukb22828_c${chr}_b0_v3.bgen" \
  -iin="${IMP_PREFIX}/ukb22828_c${chr}_b0_v3.sample" \
  -iin="${STEP1_OUT_KEEP}" \
  -iin="${SNP_QC_STEP2}" \
  -iin="${PHENO}" \
  -iin="${COV}" \
  -iin="${STEP1_OUT}fis_step_1_pred.list" \
  -icmd="bash -lc 'set -euo pipefail
  
  mkdir -p /home/dnanexus/out/out;${DOWNLOAD} 
  
  plink2 --bgen ukb22828_c${chr}_b0_v3.bgen ref-first \
  --sample ukb22828_c${chr}_b0_v3.sample \
  --extract ukb_imp.lHRC.v2.snpist \
  --keep step1_keep.txt \
  --make-pgen \
  --out ukb22828_c${chr}.HRC

  regenie --step 2 \
  --pgen ukb22828_c${chr}.HRC \
  --phenoFile fis.all.complete.pheno \
  --covarFile fis.cov \
  --pred fis_step_1_pred.list \
  --phenoColList imputed.fis2,fis2,all.fis2,iq2,imputed.iq2,all.iq2,earliest.iq,imputed.earliest.iq,all.earliest.iq2,mean.iq,imputed.mean.iq,all.mean.iq2,mean.iq.scaled,MEGAV1,MEGAV2,impvA,impvB,all.iq2.never,all.iq2.ever \
  --interaction E1 \
  --rare-mac 0 \
  --bsize 200 \
  --threads 16 \
  --gz \
  --out fis_step_2_chr${chr}

  rm ukb22828*
  rm *.loco.gz
  
  '" \
  --name "Regenie Step 2: HRC chr${chr}" \
  --instance-type "${INSTANCE}" \
  --destination "${STEP2_OUT}" \
  --priority high \
  --brief -y)
done

dx wait "${STEP2}"
echo "Step 2 complete."