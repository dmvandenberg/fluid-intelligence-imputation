############################################################### 
# Perform snipar analysis 
# 
# Run this script on Snellius/your local machine # requires: 
# -- dx toolkit on local machine 
# -- snipar-env.tar.gz on RAP 
# -- input pheno file, genotype data and covariates 
# 
###############################################################

# Snellius specific required modules
module load 2025
module load Python/3.13.1-GCCcore-14.2.0
module load SciPy-bundle/2025.06-gfbf-2025a

dx login

# IMPORTANT: project *name* or ID + trailing colon
project="UKB_IQ:"           # adjust if your project is named differently

# STATIC INPUT
SNIPAR="${project}/software/snipar/snipar-env-patch.tar.gz"
SNIPAR_REF="${project}/software/snipar"
KING="${project}/software/king"

IMP_GENO_FOLDER="${project}/Bulk/Imputation/UKB imputation from genotype/"
IMP_GENO_PREFIX="ukb22828_c"
OBS_GENO_FOLDER="${project}/Bulk/Genotype Results/Genotype calls/"
OBS_GENO_PREFIX="ukb22418_c"

# USER INPUT
AGESEX="${project}/software/snipar/agesex.txt" # file with EID, FID, age, sex
QC_SNP="${project}/reference/SNP-lists/HM3_QC_INFO98.snplist"
PHENO="${project}/input/fis.all.complete.pheno"
COV="${project}/input/fis.cov"
RAP_GWAS_DEST="${project}/output/fgwas-results"

# ==== Run King to identify first degree relatives ====

IN_BED=()
IN_BIM=()
IN_FAM=()

for chr in {1..22}; do
  IN_BED+=("-iin=${OBS_GENO_FOLDER}${OBS_GENO_PREFIX}${chr}_b0_v2.bed")
  IN_BIM+=("-iin=${OBS_GENO_FOLDER}${OBS_GENO_PREFIX}${chr}_b0_v2.bim")
  IN_FAM+=("-iin=${OBS_GENO_FOLDER}${OBS_GENO_PREFIX}${chr}_b0_v2.fam")
done

# SAK command steps:
# 1. make merge list to merge bed files into 1
# 2. Prune for LD 200 50 0.1
# 3. Export bed for pruned set of SNPs
# 4. Run King to identify 1st degree relatives
dx run swiss-army-knife \
    "${IN_BED[@]}" \
    "${IN_BIM[@]}" \
    "${IN_FAM[@]}" \
    -iin="${KING}" \
    -icmd="bash -lc '
    set -euo pipefail
    cp king   /tmp/king;   chmod +x /tmp/king

    # Merge separe chromosomes
    : > merge_list.txt
    for chr in {2..22}; do
        echo \"${OBS_GENO_PREFIX}\${chr}_b0_v2\" >> merge_list.txt
    done

    plink \
    --bfile ${OBS_GENO_PREFIX}1_b0_v2 \
    --merge-list merge_list.txt \
    --make-bed \
    --out ${OBS_GENO_PREFIX}_merged \
    --threads 16

    # LD prune (tune r2 threshold if you want; 0.1 is a common starting point)
    plink2 --bfile ${OBS_GENO_PREFIX}_merged \
      --indep-pairwise 200 50 0.1 \
      --out pruned \
      --threads 16

    # Create pruned BED
    plink2 --bfile ${OBS_GENO_PREFIX}_merged \
      --extract pruned.prune.in \
      --make-bed \
      --out ${OBS_GENO_PREFIX}_merged_pruned \
      --threads 16

    # KING on pruned set
    /tmp/king -b ${OBS_GENO_PREFIX}_merged_pruned.bed \
      --bim ${OBS_GENO_PREFIX}_merged_pruned.bim \
      --fam ${OBS_GENO_PREFIX}_merged_pruned.fam \
      --related --degree 1 \
      --prefix king_deg1 \
      --cpus 16
      
    rm ${OBS_GENO_PREFIX}*
    rm pruned.prune.in
    rm pruned.prune.out'" \
    --instance-type mem2_ssd1_v2_x16 \
    --name="King: 1st degree" \
    --destination "${SNIPAR_REF}/" \
    --priority high \
    --yes

# ==== Run IBD.py to compute IBD segments and LD scores ==== #
# Run only on observed genotypes first for IBD segments
# Steps:
# 1. Unpack snipar env
# 2. Extract HM3 SNPs
# 3. Run IBD.py
#     -- observed for IBD files
#     -- imputed for LD files
for chr in {1..22}; do
  dx run swiss-army-knife \
    -iin="${SNIPAR}" \
    -iin="${OBS_GENO_FOLDER}${OBS_GENO_PREFIX}${chr}_b0_v2.bed" \
    -iin="${OBS_GENO_FOLDER}${OBS_GENO_PREFIX}${chr}_b0_v2.bim" \
    -iin="${OBS_GENO_FOLDER}${OBS_GENO_PREFIX}${chr}_b0_v2.fam" \
    -iin="${SNIPAR_REF}/sexage.txt" \
    -iin="${COV}" \
    -iin="${SNIPAR_REF}/king_deg1.kin0" \
    -iin="${QC_SNP}" \
    -icmd="bash -lc '
      set -euo pipefail

      # 1) Unpack snipar env
      tar -xzf snipar-env-patch.tar.gz
      cd snipar-env-patch
      ./bin/conda-unpack || true
      cd ..

      # 2) Make output dir
      mkdir -p ibd_out_obs
      
      plink2 --bfile ${OBS_GENO_PREFIX}${chr}_b0_v2 \
      --extract HM3_QC_INFO98.snplist \
      --keep fis.cov \
      --make-bed \
      --out ${OBS_GENO_PREFIX}${chr}_b0_v3_hm3_0.98

      export NUMBA_DISABLE_JIT=1
      ./snipar-env-patch/bin/python ./snipar-env-patch/bin/ibd.py \
        --bed ${OBS_GENO_PREFIX}${chr}_b0_v3_hm3_0.98 \
        --king king_deg1.kin0 \
        --agesex sexage.txt \
        --chrom ${chr} \
        --out ibd_out_obs/

      rm -rf snipar-env-patch
      rm snipar-env-patch.tar.gz
      rm king_deg1.kin0
      rm sexage.txt
      rm ${OBS_GENO_PREFIX}*
      rm HM3_QC_INFO98.snplist
      rm *.cov
    ' " \
    --instance-type mem3_ssd1_v2_x2 \
    --name="Snipar IBD: chr${chr} observed" \
    --destination "${SNIPAR_REF}/" \
    --priority high \
    --yes
done

# expensive, avoid running unless crucial
for chr in {1..22}; do
  # chr 1-12 mem3_ssd2_v2x8
  # chr 13-22 mem3_ssd2_v2x4
  if (( chr >= 1 && chr <= 12 )); then
    MEM="mem3_ssd2_v2_x8"
  elif (( chr >= 13 && chr <= 22 )); then
    MEM="mem3_ssd2_v2_x4"
  fi

  dx run swiss-army-knife \
      -iin="${SNIPAR}" \
      -iin="${IMP_GENO_FOLDER}${IMP_GENO_PREFIX}${chr}_b0_v3.bgen" \
      -iin="${IMP_GENO_FOLDER}${IMP_GENO_PREFIX}${chr}_b0_v3.sample" \
      -iin="${SNIPAR_REF}/sexage.txt" \
      -iin="${COV}" \
      -iin="${SNIPAR_REF}/king_deg1.kin0" \
      -iin="${QC_SNP}" \
      -icmd="bash -lc '
        set -euo pipefail

        # 1) Unpack snipar env
        tar -xzf snipar-env-patch.tar.gz
        cd snipar-env-patch
        ./bin/conda-unpack || true
        cd ..

        # 2) Make output dir
        mkdir -p ibd_out_imp
        
        plink2 --bgen ${IMP_GENO_PREFIX}${chr}_b0_v3.bgen ref-first \
        --sample ${IMP_GENO_PREFIX}${chr}_b0_v3.sample \
        --extract HM3_QC_INFO98.snplist \
        --keep fis.cov \
        --make-bed \
        --out ${IMP_GENO_PREFIX}${chr}_b0_v3_hm3_0.98

        export NUMBA_DISABLE_JIT=1
        ./snipar-env-patch/bin/python ./snipar-env-patch/bin/ibd.py \
          --bed ${IMP_GENO_PREFIX}${chr}_b0_v3_hm3_0.98 \
          --king king_deg1.kin0 \
          --agesex sexage.txt \
          --chrom ${chr} \
          --ld_out \
          --out ibd_out_imp/

        rm -rf snipar-env-patch
        rm snipar-env-patch.tar.gz
        rm king_deg1.kin0
        rm sexage.txt
        rm ${IMP_GENO_PREFIX}*
        rm HM3_QC_INFO98.snplist
        rm *.cov
      ' " \
      --instance-type "${MEM}" \
      --name="Snipar IBD: chr${chr} imputed" \
      --destination "${SNIPAR_REF}/" \
      --priority high \
      --yes
done

# ==== Run impute.py to impute genotypes to form trio's ==== #
# Run only on observed genotypes IBD segments
# Steps:
# 1. Unpack snipar env
# 2. Extract HM3 SNPs from imputed genotypes
# 3. Run impute.py
for chr in {1..22}; do
  dx run swiss-army-knife \
    -iin="${SNIPAR}" \
    -iin="${IMP_GENO_FOLDER}${IMP_GENO_PREFIX}${chr}_b0_v3.bgen" \
    -iin="${IMP_GENO_FOLDER}${IMP_GENO_PREFIX}${chr}_b0_v3.sample" \
    -iin="${SNIPAR_REF}/sexage.txt" \
    -iin="${COV}" \
    -iin="${SNIPAR_REF}/king_deg1.kin0" \
    -iin="${QC_SNP}" \
    -iin="${SNIPAR_REF}/ibd_out_obs/_chr_${chr}.ibd.segments.gz" \
    -icmd="bash -lc '
      set -euo pipefail

      # 1) Unpack snipar env
      tar -xzf snipar-env-patch.tar.gz
      cd snipar-env-patch
      ./bin/conda-unpack || true
      cd ..

      # 2) Make output dir
      mkdir -p HM3_098
      mkdir -p HDF

      plink2 --bgen ${IMP_GENO_PREFIX}${chr}_b0_v3.bgen ref-first \
        --sample ${IMP_GENO_PREFIX}${chr}_b0_v3.sample \
        --chr ${chr} \
        --extract HM3_QC_INFO98.snplist \
        --keep fis.cov \
        --make-bed \
        --out HM3_098/${IMP_GENO_PREFIX}${chr}_b0_v3_hm3_0.98.EU \
        --threads 16 \
      
      # 3) Impute parental genotypes
      export NUMBA_DISABLE_JIT=1
      ./snipar-env-patch/bin/python ./snipar-env-patch/bin/impute.py \
        --bed HM3_098/${IMP_GENO_PREFIX}${chr}_b0_v3_hm3_0.98.EU \
        --ibd _chr_${chr}.ibd \
        --king king_deg1.kin0 \
        --agesex sexage.txt \
        --threads 16 \
        --out HDF/chr${chr}

      rm -rf snipar-env-patch
      rm snipar-env-patch.tar.gz
      rm king_deg1.kin0
      rm HM3_QC_INFO98.snplist
      rm sexage.txt
      rm ${IMP_GENO_PREFIX}${chr}*
      rm 25PCs.age.age2.cov
      rm *.ibd*
      rm -rf HM3_098

    ' " \
    --instance-type mem3_ssd1_v2_x16 \
    --name="Snipar impute: chr${chr}" \
    --destination "${SNIPAR_REF}/" \
    --priority high \
    --yes
done

# ==== Run gwas.py to perform within-family GWAS ==== #
# Uses .hdf files from impute.py and imputed genotypes as input
# Steps:
# 1. Unpack snipar env
# 2. Extract HM3 SNPs from imputed genotypes
# 3. Run gwas.py for desired phenotypes in phenofile

# change according to pheno names in dedicated file
for pheno_name in  "fis2" "all.iq2" "earliest.iq" "all.earliest.iq2" "mean.iq" "all.mean.iq2" "MEGAV1" "MEGAV2"; do
  for chr in {1..22}; do
    dx run swiss-army-knife \
      -iin="${SNIPAR}" \
      -iin="${PHENO}" \
      -iin="${COV}" \
      -iin="${SNIPAR_REF}/HDF/chr${chr}.hdf5" \
      -iin="${QC_SNP}" \
      -iin="${IMP_GENO_FOLDER}${IMP_GENO_PREFIX}${chr}_b0_v3.bgen" \
      -iin="${IMP_GENO_FOLDER}${IMP_GENO_PREFIX}${chr}_b0_v3.sample" \
      -icmd="bash -lc '
        set -euo pipefail

        # 1) Unpack snipar env
        tar -xzf snipar-env-patch.tar.gz
        cd snipar-env-patch
        ./bin/conda-unpack || true
        cd ..

        # 2) Make output dir
        mkdir -p ${pheno_name}

        # 3) subset genotypes
        plink2 --bgen ${IMP_GENO_PREFIX}${chr}_b0_v3.bgen ref-first \
        --sample ${IMP_GENO_PREFIX}${chr}_b0_v3.sample \
        --chr ${chr} \
        --extract HM3_QC_INFO98.snplist \
        --keep fis.cov \
        --make-bed \
        --out ${IMP_GENO_PREFIX}${chr}_b0_v3_hm3_0.98.EU \
        --threads 16 \

        # 5) Run SNIPAR gwas
        export NUMBA_DISABLE_JIT=1
        ./snipar-env-patch/bin/python ./snipar-env-patch/bin/gwas.py \
          fis.all.complete.pheno \
          --bed ${IMP_GENO_PREFIX}${chr}_b0_v3_hm3_0.98.EU \
          --imp chr${chr} \
          --out ${pheno_name}/chr${chr} \
          --phen ${pheno_name} \
          --cov fis.cov \
          --no_hdf5_out \

        rm -rf snipar-env-patch
        rm snipar-env-patch.tar.gz
        rm fis.cov
        rm chr${chr}.hdf5
        rm ${IMP_GENO_PREFIX}${chr}*
        rm HM3_QC_INFO98.snplist
        rm fis.all.complete.pheno
      ' " \
      --instance-type mem3_ssd1_v2_x16 \
      --name="Snipar GWAS: ${pheno_name}, chr ${chr}" \
      --destination "${RAP_DEST}" \
      --priority high \
      --yes
  done
done

# Download LD scores and sumstats and run correlate.py through parse-snipar.sh
