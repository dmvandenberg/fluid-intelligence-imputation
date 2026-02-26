# Fluid intelligence imputation in UK Biobank
In this repository we provide a general lay summary and FAQ along with scripts and GWAS summary statistics used in the study 'Imputation of fluid intelligence scores reduces ascertainment bias and increases power for analyses of common and rare variants' ([preprint](https://www.medrxiv.org/content/10.1101/2025.06.18.25329418v1)). The scripts cover phenotype preparation, SoftImpute-based imputation, GWAS, and downstream analyses. 

## General Lay Summary
Research on the genetics of intelligence offers both scientific opportunities and ethical responsibilities. We carefully considered the goal of this research, which was, broadly, to better understand the genetics of measures of intelligence in order to eventually inform our understanding of the causes of clinically recognised intellectual disability and developmental delay. Our research was approved by UK Biobank. In addition, we have prepared this lay summary and FAQ to help broader audiences engage with the motivation for this research and the results.  

This study explores how we can improve genetic research on cognitive ability when a large part of the data is missing. Due to its large sample size, most of the statistical power in genetic research on cognitive ability so far has come from UK Biobank. About 40 percent of UK Biobank participants, however, never completed the fluid intelligence test, which limits what previous studies could discover. The missing data also introduces bias, since people who took the test tend to be more educated.

To address this, we estimated the missing fluid intelligence test scores using information from a wide range of related traits, including physical and mental health, behavior, and socio-economic outcomes. We then combined these estimated scores with the actual test scores to conduct the largest genetic analyses of cognitive ability to date. The analyses were restricted to UK Biobank participants of genetically inferred European ancestry.

This approach led to the discovery of many more genetic variants associated with cognition, including rare variants in genes linked to neurodevelopment. These rare variant associations would not have been detectable without the increased sample size provided by imputation, and many are located in genes known to affect brain development. Importantly, the estimated fluid intelligence scores closely matched the genetic patterns of the real test scores and helped reduce bias caused by missing data.

## Frequently Asked Questions (FAQ)
###Why study the genetics of intelligence?
Firstly, it is important to note that we can only study intelligence through measurable and imperfect proxies, such as specific cognitive tests, and not as some idealized or fixed concept. When we talk about studying the “genetics of intelligence/cognitive ability”, we really mean we are studying the genetics of specific measures of intelligence/cognitive ability.

Understanding the genetic and environmental factors that shape cognitive differences can help us better understand certain health conditions, as well as how cognition relates to social and health inequalities. Many of the same genetic variants that influence cognitive ability also affect the risk for neurodevelopmental conditions. Socio-economic outcomes like educational attainment (how long do people spend in education) share roughly half of their genetic influences with intelligence. Because cognitive ability plays such a central role in many modern societies, its relationship with genetics is sometimes controversial and vulnerable to misinterpretation. We believe this is important research that should be carried out carefully, transparently, and ethically.

### What did the study find?
- We increased the number of people in our genetic analyses from about 270,000 to over 450,000.
- This led to the discovery of more common genetic variants linked to cognitive ability, rising from 390 to 550.
- We detected 26 genes where rare damaging variants are linked to lower cognitive ability, including 14 that were already known to play a role in neurodevelopmental conditions.
- Most people carrying these rare variants in UK Biobank do not have a diagnosis of a neurodevelopmental disorder. This suggests these rare genetic effects can influence cognition even without causing people to reach a clinical threshold for diagnosis.
- The improved data led to more powerful genetic predictors of cognitive ability in other studies.


### What does “imputing intelligence” mean?
We used a statistical model that learns how a wide range of physical, psychological, and socio-economic outcomes relate to intelligence test scores in people who completed the test. We then used that model to estimate likely scores for individuals who did not take it. This approach is commonly used in science to handle missing data. It helps increase statistical power and reduce bias, especially when participation in the test is not random.

### Are there sensitivities surrounding the imputation of cognitive scores for participants who chose not to complete the test?
We acknowledge that imputing cognitive scores for missing participants touches on sensitive ground. In UK Biobank, approximately 60% of participants completed the fluid intelligence test at least once, and non-participation arose from various factors, including study design (e.g. the test was only introduced into the ‘baseline assessment’ part way through recruitment), optional participation, and time constraints. Importantly, participants who formally withdraw their consent are excluded from the dataset, and we regularly updated our data to reflect these exclusions. Imputation is a standard scientific method used to handle missing data and increase statistical power. This approach not only boosts power but also helps reduce biases that arise when analyses include only participants with complete data, who often differ systematically from those with missing data. Other approaches are routinely applied to infer characteristics not directly reported by participants, such as genetic ancestry or environmental exposures. Importantly, our imputed fluid intelligence scores are used solely for research purposes at the group level and are not intended to make predictions about any individual participant.

### Does this study mean intelligence is genetic?
It was well-established before this study that genes influence cognitive ability, but also that they are only part of the story. The environment, including education, parenting, societal context and health-related factors, also plays a major role. Genetic effects can interact with and correlate with environmental factors. This study helps identify some of the genes that may be involved, but it does not mean that intelligence is fixed or determined at birth.

### How generalizable is this research?
This study was conducted with individuals of genetically-inferred European ancestries only. That is because most of the available genetic data comes from this group. Results from one  group of people are difficult to apply directly to others, not only because of cultural and environmental differences, but also for technical reasons. The patterns of correlation between genetic variants can differ across ancestries, which can affect statistical results even when the underlying genetic effects are similar. This is a serious limitation in many genetic studies, and the field needs to become more inclusive by collecting and analyzing data from more diverse populations.

### Can this research be used to study population differences in measures of intelligence?
While this question lies well beyond the scope of our study, we address it here because it sits at the heart of some of the public controversy around genetic research on cognitive traits. We are aware that fringe researchers have attempted to (mis)use genetic results and summary statistics in this area, and we believe it is important to clarify why such comparisons are scientifically problematic and should be approached with caution.

Firstly, human ancestries are complex, shaped by overlapping migrations, histories, and admixture. Categories like “European”, “South Asian”, or “Black African” are socially and historically constructed and do not always map cleanly onto underlying genetic patterns. This study only includes individuals of European ancestries and cannot be used to make claims about differences between ancestry groups. Applying results across populations is limited by both social and technical factors.

Even if average polygenic indices differ between populations, that does not mean there are true genetic differences in the traits themselves. One reason these indices might differ is due to differences in the pattern of correlations between genetic variants across ancestries; a difference in polygenic indices between ancestries does not necessarily imply systematic differences in the frequency of genetic variants that have causal effects on the trait in question (i.e. the difference in polygenic indices can exist even if the underlying biology of the trait is the same between ancestries). Furthermore, there is evidence that genetic effects on cognitive traits are sensitive to social and economic conditions, which differ between groups and are shaped by history, inequality, and discrimination. Finally, natural selection can complicate comparisons. For traits that are kept close to an optimal value by natural selection (e.g., through stabilizing selection), different combinations of genetic variants can evolve differently in different populations to maintain the same average trait level. This means that polygenic indices based on one population might suggest differences that are not actually there or miss real differences that do exist.

Because of these issues, genetic comparisons of cognitive traits between populations are both very complex to execute robustly and highly sensitive to misinterpretation. 

### Are these findings relevant for clinical or educational use?
These findings are meant to improve our scientific understanding of cognition, not to make individual predictions. Although it is technically possible to make predictions of measures of intelligence using measures like polygenic indices, we do not use polygenic indices to assess individuals in this paper and do not recommend doing so. Polygenic indices based on this research are noisy and reflect subtle differences that are only useful when studying large groups. They remain inferior to direct measurements of cognition-related traits, such as standardized exam results in school systems. They are neither accurate nor reliable enough for use in personal decision-making in health, education, or employment.

The same caution applies to the rare genetic variants identified in this study. While some of the genes that show rare variant associations with fluid intelligence are known to play a role in neurodevelopmental conditions, many people who carry rare variants in these genes do not have a recorded clinical diagnosis in UK Biobank. Our findings highlight potential difficulties in interpreting rare variants in these genes found in patients with a rare neurodevelopmental condition, particularly if those variants have been inherited from clinically unaffected family members. This problem of interpreting variants like this with “incomplete penetrance” is likely to become an increasing challenge for clinical genetics over the coming years. The penetrance of these variants for neurodevelopmental conditions needs to be evaluated in very large, minimally ascertained cohorts (e.g. population-wide birth cohorts) to inform their interpretation in clinical contexts. 

### Can this research be used for embryo selection?
Many countries, including the UK, have strict regulations on the use of polygenic indices in embryo selection. However, some countries allow their use under certain circumstances.  This research was not designed for embryo selection, and applying it in that context raises scientific and ethical concerns. Polygenic indices have limited predictive power within families (i.e. the scenario in which they would be used when selecting one of a set of embryos produced from the same couple), especially for traits like intelligence. Predictive power may increase as discovery sample sizes grow, particularly if missing cognitive phenotypes can be recovered using similar imputation approaches across multiple biobanks. If predictive power does improve, it would increase the stakes of trade-offs and unintended effects on other traits. Even when predictive power is modest, embryo selection could potentially still shift correlated risks in unintended directions. For example, the genetic variants that are associated with higher intelligence also tend to increase the risk for autism. More importantly, these indices reflect a mix of biological and social influences, many of which are shaped by the current environment and may change over time. Because complex traits are influenced by many overlapping biological and environmental pathways, the long-term effects of selecting embryos in this way are unpredictable.

### What is the next step?
We hope other researchers will use our methods to improve the quality of genetic research in large datasets. We also intend to use our results to see if they help us discover more genes associated with neurodevelopmental conditions in clinical cohorts. Future studies should apply these approaches in more diverse populations, and continue to explore how genes and environments interact to shape cognitive outcomes and neurodevelopmental conditions.

## Computation details
Here we list scripts used in the paper "Imputation of fluid intelligence scores reduces ascertainment bias and increases power for analyses of common and rare variants"
for a detailed description of analyses and decisions please refer to our paper.

### Imputation
Scripts regarding imputation

- **[prepare_FI_pheno.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/prepare_FI_pheno.R)**: is the script used to standardize FI measured and configure covariate files for GWAS with the appropriate values (age at time of measurement). This script requires UKB variables 34-0.0 and 52-0.0, and all instances of 20016, 20131, 21003, 20135. In addition it requires PC loadings for 25 ancestry PCs, and the sex and array used for genotyping for all individuals.
- **[FI_imputation.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/FI_imputation.R)**: is used for imputation. It leverages softImpute and requires access to various UKB datafield to impute FI. We have evaluated 2 sets of imputation variables for which the UKB variable ID's can be found [here](https://github.com/dmvandenberg/UKB-FI-Imputation/tree/main/imputation/imputation-variables).
- **[GWAS-by-subtraction](https://github.com/dmvandenberg/UKB-FI-Imputation/tree/main/imputation/gwas-by-subtraction)**: Consists of 3 main scripts: 1) [Prepare_IQ_ImpIQ.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/gwas-by-subtraction/Prepare_IQ_ImpIQ.R), applies GSEM munge and does some cleaning of SNPs; 2) [Subtract_IQ_ImpIQ.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/gwas-by-subtraction/Subtract_IQ_ImpIQ.R), performs GWAS by subtraction, parallelizing chunks. 3) [extract_gsem_blocks.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/gwas-by-subtraction/extract_gsem_blocks.R), extract output and compile GWAS.
- **[rescale_combine.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/rescale_combine.R)**: is used to mega-analyze measured and imputed FI values after standardizing them separately.
- **[combine.metal](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/imputation/combine.metal)**: is used to meta-analyze GWAS of measured and imputed FI.

### Genetic analysis

- **[gwas.job](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/genetic-analyses/gwas.job)**: is the script used to run a GWAS using [GCTA](https://yanglab.westlake.edu.cn/software/gcta/). This requires access to genotype data, a sparse GRM, and one of the .pheno files created by our imputation scripts along with accompanying .qcov file
- **[COGENT_meta.metal](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/genetic-analyses/COGENT_meta.metal)**: is the [metal](https://genome.sph.umich.edu/wiki/METAL_Documentation) meta-analysis script used to meta-analyze the GWAS of COGENT and Combined FIS
  
#### Common variant analyes

Genetic correlations and SNP heritabilities were computed using [LDSC](https://github.com/bulik/ldsc). The scripts we used for this are:

- **[munge.job](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/genetic-analyses/munge.job)**: Script to format following LDSC standards
- **[h2_rg.job](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/genetic-analyses/h2_rg.job)**: To compute SNP heritability and genetic correlations with reference GWAS

#### Rare variant analyses

1. Quality control of rare genetic variants
   We used applets designed by Eugene Gardner to process the quality control of rare variants. For details, please refer to https://github.com/mrcepid-rap. specifically, we used mrc-splitbcf, mrc-filterbcf, mrc-makebgen and mrc-collapsevariants to obtain the variant   files for regenie. We filtered the variants by minor allele frequency (MAF), retaining only those with a MAF lower than 0.001%. We also annotated two additional variant effect predictor MPC and alphamissense to the missense variants using in-house written script [annotate-score.sh](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/rare-variant-analyses/annotate-score.sh).

   For PTVs, only high-confidence PTVs defined by LOFTEE were retained. For missense variants, a damaging missense variant set was created by including variants (Alphamissense >= 0.56, & REVEL >= 0.5 & MPC >=2).  

2. format the variant file for regenie
   We specifically use <output_prefix>..STAAR.variants_table.tsv from the Output Tarball of mrc-collapsevariants. We used [UKB_fileformat.R](https://github.com/dmvandenberg/UKB-FI-Imputation/blob/main/rare-variant-analyses/UKB_fileformat.R) to produce the setlist and annotation files required to run gene-based testing in Regenie.

3. Running regenie on UKB research analysis platform (UKB-RAP)
   After files were prepared for regenie, we used two steps to run gene-based testing in Regenie on UKB-RAP. We used swiss-army-knife to run regenie.
   a) Pre-filtering before running regenie.

   Merge common genetic variant files
   `dx run swiss-army-knife -iin="default.txt" -icmd='cp /project/Bulk/Imputation/UKB\ imputation\ from\ genotype/ukb22418_c[1-9]* . ;ls *.bgen | sed -e 's/.bgen//g'> files_to_merge.txt;plink --merge-list files_to_merge.txt --make-bed --autosome-xy --out ukb22828_c1_22_v2_merged;rm files_to_merge.txt;' --tag='Step1' --instance-type "mem1_ssd1_v2_x16" --destination='/project/Data/' --brief --yes`

   Quality control the SNPs
   `dx run swiss-army-knife -iin="/project/Data/ukb22418_c1_22_v2_merged.bed" -iin="/project/Data/ukb22418_c1_22_v2_merged.bim" -iin="/project/Data/ukb22418_c1_22_v2_merged.fam" -icmd="plink2 --bfile ukb22418_c1_22_v2_merged --autosome --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 --write-snplist --write-samples --no-id-header --out qc_pass" --tag="Array QC" --instance-type "mem2_ssd2_v2_x96" --destination="/project/Data/" --brief --yes`

   b) Example: regenie step1

   Input files:
   - genetic data file: bed/bim/fam file
   - phenotype data file: ukb_phenotype.txt including focal phenotype and covariates
   - covariate data file: same with phenotype data file
   - optional: list of variants to keep (qc_pass.snplist)
   - output file: Predictions files (.loco.gz)
  
   `dx run swiss-army-knife -iin="ukb22418_c1_22_v2_merged.bed" -iin="ukb22418_c1_22_v2_merged.bim" -iin="ukb22418_c1_22_v2_merged.fam" -iin="qc_pass.snplist" -iin="ukb_phenotype.txt" -icmd="regenie --step 1 --bed ukb22418_c1_22_v2_merged --phenoFile ukb_phenotype.txt --covarFile ukb_phenotype.txt --extract qc_pass.snplist --phenoColList fluid.intel --covarColList age,birth.year,pca{1:25},age.squared,agebysex,age.squaredbysex,sex --catCovarList center --maxCatLevels 30 --lowmem --bsize 1000 --qt --apply-rint --gz --verbose --threads 16 --out fis_step1" --instance-type "mem1_ssd1_v2_x16" --destination='/result/' --brief --yes`

    c) Example: regenie step2 (eg.chr1)

   Input files:
   - genetic data file: BGEN file and Sample file corresponding to input BGEN file
   - phenotype data file: ukb_phenotype.txt including focal phenotype and covariates
   - covariate data file: same with phenotype data file
   - prediction files: ...pred.list and ...loco.gz from step1 results
   - annotation input files: ...setlist.txt and ...annotations.txt (see above for file preparations)
   - mask files: mask.txt specifies which annotation categories should be combined into masks
   - sample files: ..BOLT.sample (BOLT-ready bgen sample file of per-gene 'genotypes' from mrc-collapsevariants
   - output file: regenie results files
  
   `dx run swiss-army-knife -iin="/filtered-bgen/chr1.filtered.bgen" -iin="/project/Data/HC-PTV-chr1.22.BOLT.sample"   -iin="/project/Data/ukb_phenotype.txt" -iin="/project/Data/fis_step1_pred.list" -iin="/project/Data/fis_step1_1.loco.gz" -iin="/project/step2_lof/chr1_setlist.txt" -iin="/project/step2_lof/chr1_annotations.txt" -iin="/project/step2_lof/mask.txt" -icmd="regenie --step 2 --bgen chr1.filtered.bgen --sample HC-PTV-chr1.22.BOLT.sample --phenoFile ukb_phenotype.txt --covarFile ukb_phenotype.txt --phenoColList fluid.intel --covarColList age,birth.year,pca{1:25},age.squared,agebysex,age.squaredbysex,sex --catCovarList center --maxCatLevels 30 --firth --approx --pred fis_step1_pred.list  --anno-file chr1_annotations.txt --set-list chr1_setlist.txt --mask-def mask.txt --aaf-bins 1 --vc-tests skato-acat,acato-full --bsize 200 --qt --apply-rint --threads 2 --gz --verbose --out chr1_ptv" --tag='Step1' --instance-type "mem1_ssd1_v2_x16" --destination='/project/Data/' --brief --yes`
