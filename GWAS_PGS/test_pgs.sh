#!/usr/bin/env bash

#SBATCH -J PGS
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=75GB
#SBATCH --partition=tier2q
#SBATCH -o /home/mcostantino/job-output/run_PGS/PGS_PGSC_%A_%a.out
#SBATCH -e /home/mcostantino/job-output/run_PGS/PGS_PGSC_%A_%a.err
#SBATCH --array=1-39

now=$(date)
# Record job start time
start_time=$(date +%s)
echo "Date run: $now"


lambda=1
GENOS="/gpfs/data/ukb-share/genotypes/pop_genos/"
PHENOS="/gpfs/data/ukb-share/extracted_phenotypes/"
LOCAL_PLINK2="/ess/home/home1/mcostantino/plink2"
TRAIN_POP="whitebrit"
TEST_POP="whitebrit"
JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results/"
pop_path="${GENOS}${TEST_POP}/"
genos="${GENOS}${TEST_POP}/"
pop_file="ukb_chr1-22_${TEST_POP}"
PGS="/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out/"
CLUMPED_GENOS="/scratch/mcostantino/pgs_output/clumped_genos/"
THRESH=(0.0000000001 0.00000001 0.000001 0.0001 0.001 0.005 0.01 0.05 0.1 0.5)
PHENOARRAY=(
  "alcohol_intake_frequency"
  "arm_fat-free_mass_left"
  "arm_fat-free_mass_right"
  "basophill_count"
  "birth_weight"
  "bmi"
  "calcium"
  "cholesterol"
  "creatinine"
  "diastolicbp_auto"
  "eosinophill_count"
  "glucose"
  "hdl"
  "hba1c"
  "height"
  "hip_circumference"
  "igf-1"
  "ldl"
  "leukocyte_count"
  "mean_corpuscular_volume"
  "platelet_count"
  "pulse_rate"
  "rbc"
  "shbg"
  "systolicbp_auto"
  "testosterone"
  "triglycerides"
  "urate"
  "urea"
  "vitamin_d"
  "waist_circumference"
  "whole_body_fat_mass"
  "whradjbmi_zhu"
  "whradjbmi_emdin"
  "fev1"
  "ea4"
  "c-reactive_protein"
  "hip_to_waist"
  "arm_fat-free_mass_avg"
)

if (( $(echo "$lambda == 0" | bc -l) )); then
	  JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results_log/"
    PGS="/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out_log/"
    CLUMPED_GENOS="/scratch/mcostantino/pgs_output/clumped_genos_log/"
fi

echo "SLURM_ARRAY_TASK_ID: $SLURM_ARRAY_TASK_ID"
echo "Calculated pheno_num: $pheno_num"
echo "PHENOARRAY length: ${#PHENOARRAY[@]}"
# Calculate the index for the selected phenotype
pheno_num=$((SLURM_ARRAY_TASK_ID-1))
selected_pheno="${PHENOARRAY[pheno_num]}"
#selected_pheno="diastolicbp_auto"

echo "$selected_pheno"
if [[ ! $selected_pheno ]]
  then
    echo "no pheno present"
    exit
fi

pgs_path="${PGS}${selected_pheno}/"

mkdir -p ${CLUMPED_GENOS}
mkdir -p ${pgs_path}
mkdir -p ${JOINED_GWAS}snp_pval/

# REMOVE THIS LOOP
echo "Building std ${selected_pheno} PGS"
# clump snps
output_file="${CLUMPED_GENOS}${pop_file}_${selected_pheno}"
gwas_file="${JOINED_GWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}.assoc.linear"	
valid_snp_file="${CLUMPED_GENOS}${pop_file}_${selected_pheno}.valid.snp"
if [[ ! -f "$output_file".clumps ]]; then
  # clump snps for prs
  ${LOCAL_PLINK2} \
    --bfile ${pop_path}${pop_file}_QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump "$gwas_file" \
    --out "$output_file"
  echo "$selected_pheno clump done"
	        
  # grab col 3 for rsids
  awk 'BEGIN{OFS=FS} NR!=1{print $3}' ${output_file}.clumps > "$valid_snp_file"
else
  echo "Output file already exists at $output_file"
  awk 'BEGIN{OFS=FS} NR!=1{print $3}' ${output_file}.clumps > "$valid_snp_file"
fi
  
# build pgs
pgs_out_file="${pgs_path}${selected_pheno}_pgs"
last_thresh="${THRESH[-1]}" #grab ex. thresh
# for the --score flag : 3 5 12 grabs snp ID, alt allele, & effect size estimate
# for --q-score-range : cols 3 & 15 of "$gwas_file" grabs snp IDs & pval cols
if [[ ! -f "$pgs_out_file"."$last_thresh".sscore ]]; then
  ${LOCAL_PLINK2} \
    --bfile ${pop_path}${pop_file}_QC \
    --keep ${genos}${pop_file}_QC_test \
    --score "$gwas_file" 3 5 12 header cols=+scoresums \
    --q-score-range std_range_list "$gwas_file" 3 15 header \
    --extract "$valid_snp_file" \
    --out "$pgs_out_file"
  echo "Std PGS complete"
  #add back to cut pop to 1K
  #--keep ${pop_path}/${pop_file}_QC_1000.fam \
else
  echo "std pgs already exists"
fi
		
echo "test PGS done"

# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds" 