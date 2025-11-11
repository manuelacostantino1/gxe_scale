#!/usr/bin/env bash

#SBATCH -J ldsc
#SBATCH --time=20:00:00
#SBATCH --mem=200GB
#SBATCH --partition=tier2q
#SBATCH -o /home/mcostantino/job-output/ldsc/ldsc.out
#SBATCH -e /home/mcostantino/job-output/ldsc/ldsc.err

# I think these two commands transform the summary statistics into the correct formats
pheno="testosterone"
LSDC_PATH="/gpfs/data/ukb-share/dahl/manuela/scale_project/descaler/ldsc"
INTERMEDIATE="/scratch/mcostantino/descaler/intermediate_files"
LD_DIR="/scratch/mcostantino/descaler/intermediate_files/LDscore"
#downloaded the ld references from https://zenodo.org/records/10515792/files/1000G_Phase3_ldscores.tgz?download=1
#downloaded hapmap3 with wget https://raw.githubusercontent.com/bulik/ldsc/master/w_hm3.snplist

module load gcc/12.1.0 miniconda3/24.9.2 ensembl-vep/114conda
cd $LSDC_PATH
source activate ldsc

# Transform GWAS output for default scale
if [ ! -f $INTERMEDIATE/default_${pheno}.sumstats.gz ]; then
  $LSDC_PATH/munge_sumstats.py \
    --sumstats /gpfs/data/ukb-share/dahl/manuela/joined_gwas_results/${pheno}/chr1-22_whitebrit_${pheno}.assoc.linear \
    --snp ID \
    --a1 A1 \
    --a2 OMITTED \
    --p P \
    --signed-sumstats BETA,0 \
    --N-col OBS_CT \
    --out $INTERMEDIATE/default_${pheno} \
    --chunksize 500000 \
    --merge-alleles $INTERMEDIATE/w_hm3.snplist
fi

if [ ! -f $INTERMEDIATE/log_${pheno}.sumstats.gz ]; then
  $LSDC_PATH/munge_sumstats.py \
    --sumstats /gpfs/data/ukb-share/dahl/manuela/joined_gwas_results_log/${pheno}/chr1-22_whitebrit_${pheno}.assoc.linear \
    --snp ID \
    --a1 A1 \
    --a2 OMITTED \
    --p P \
    --signed-sumstats BETA,0 \
    --N-col OBS_CT \
    --out $INTERMEDIATE/log_${pheno} \
    --chunksize 500000 \
    --merge-alleles $INTERMEDIATE/w_hm3.snplist
fi

if [ ! -f $INTERMEDIATE/males_${pheno}.sumstats.gz ]; then
  $LSDC_PATH/munge_sumstats.py \
    --sumstats /gpfs/data/ukb-share/dahl/manuela/joined_gwas_male/${pheno}/chr1-22_whitebrit_${pheno}_males.glm.linear \
    --snp ID \
    --a1 A1 \
    --a2 OMITTED \
    --p P \
    --signed-sumstats BETA,0 \
    --N-col OBS_CT \
    --out $INTERMEDIATE/males_${pheno} \
    --chunksize 500000 \
    --merge-alleles $INTERMEDIATE/w_hm3.snplist
fi

if [ ! -f $INTERMEDIATE/females_${pheno}.sumstats.gz ]; then
  $LSDC_PATH/munge_sumstats.py \
    --sumstats /gpfs/data/ukb-share/dahl/manuela/joined_gwas_female/${pheno}/chr1-22_whitebrit_${pheno}_females.glm.linear \
    --snp ID \
    --a1 A1 \
    --a2 OMITTED \
    --p P \
    --signed-sumstats BETA,0 \
    --N-col OBS_CT \
    --out $INTERMEDIATE/females_${pheno} \
    --chunksize 500000 \
    --merge-alleles $INTERMEDIATE/w_hm3.snplist
fi

if [ ! -f $INTERMEDIATE/logfemales_${pheno}.sumstats.gz ]; then
  $LSDC_PATH/munge_sumstats.py \
    --sumstats /gpfs/data/ukb-share/dahl/manuela/joined_gwas_female_log/${pheno}/chr1-22_whitebrit_${pheno}_females.glm.linear \
    --snp ID \
    --a1 A1 \
    --a2 OMITTED \
    --p P \
    --signed-sumstats BETA,0 \
    --N-col OBS_CT \
    --out $INTERMEDIATE/logfemales_${pheno} \
    --chunksize 500000 \
    --merge-alleles $INTERMEDIATE/w_hm3.snplist
fi

if [ ! -f $INTERMEDIATE/logmales_${pheno}.sumstats.gz ]; then
  $LSDC_PATH/munge_sumstats.py \
    --sumstats /gpfs/data/ukb-share/dahl/manuela/joined_gwas_male_log/${pheno}/chr1-22_whitebrit_${pheno}_males.glm.linear \
    --snp ID \
    --a1 A1 \
    --a2 OMITTED \
    --p P \
    --signed-sumstats BETA,0 \
    --N-col OBS_CT \
    --out $INTERMEDIATE/logmales_${pheno} \
    --chunksize 500000 \
    --merge-alleles $INTERMEDIATE/w_hm3.snplist
fi

# Note: if i give it three files it will compute the rg for the first file with all subsequent files
# This command runs the LD Score regression analysis

methods=("default" "log" "males" "females" "logmales" "logfemales")

mkdir -p /gpfs/data/ukb-share/dahl/manuela/ldsc_results/${pheno}

for method1 in "${methods[@]}"; do
  for method2 in "${methods[@]}"; do
    if [[ "$method1" != "$method2" ]]; then
      ./ldsc.py \
        --rg $INTERMEDIATE/${method1}_${pheno}.sumstats.gz,$INTERMEDIATE/${method2}_${pheno}.sumstats.gz \
        --ref-ld-chr ${LD_DIR}/LDscore. \
        --w-ld-chr ${LD_DIR}/LDscore. \
        --print-snps $INTERMEDIATE/w_hm3.snplist \
        --out /gpfs/data/ukb-share/dahl/manuela/ldsc_results/${pheno}/${method1}_${method2}_${pheno}
    fi
  done
done

