#!/usr/bin/env bash

#SBATCH -J val_PGS_PGSC
#SBATCH --time=09:00:00
#SBATCH --mem-per-cpu=10GB
#SBATCH --partition=tier1q
#SBATCH -o /home/mcostantino/job-output/run_val_PGS/val_PGS_PGSC_%A_%a.out
#SBATCH -e /home/mcostantino/job-output/run_val_PGS/val_PGS_PGSC_%A_%a.err
#SBATCH --array=1-39

now=$(date)
# Record job start time
start_time=$(date +%s)
echo "Date run: $now"

lambda=0
GENOS="/gpfs/data/ukb-share/genotypes/pop_genos/"
LOCAL_PLINK2="/ess/home/home1/mcostantino/plink2"
TRAIN_POP="whitebrit"
TEST_POP="whitebrit"
VALID_POPS=("white_euro" "afr" "asn")
pop_path="${GENOS}${TEST_POP}/"
genos="${GENOS}${TEST_POP}/"
pop_file="ukb_chr1-22_${TEST_POP}"
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
    R2="/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_log_log_prediction/"
elif (( $(echo "$lambda == 1" | bc -l) )); then
    JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results/"
    PGS="/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/pgs_out/"
    R2="/gpfs/data/ukb-share/dahl/manuela/pgs_outputs/r2_out_default_log_prediction/" ## Taking in different inputs
    CLUMPED_GENOS="/scratch/mcostantino/pgs_output/clumped_genos/"
fi

best_thresh="${PGS}best_thresh_log/"
pgs_path="${PGS}valid_log/" ## Sending the output files to differnt folder than before

mkdir -p ${best_thresh}
mkdir -p ${pgs_path}

# Calculate the index for the selected phenotype
pheno_num=$((SLURM_ARRAY_TASK_ID-1))
selected_pheno="${PHENOARRAY[pheno_num]}"

for valid_pop in ${VALID_POPS[@]}; do

  pop_path="${GENOS}${valid_pop}/"
  genos="${GENOS}${valid_pop}/"
  pop_file="ukb_chr1-22_${valid_pop}"


  echo "Building std ${selected_pheno} PGSs, for ${valid_pop}"

    # grab best thresh from test pgs and use it below instead of std_range_list
    processed_pgs="${R2}processed_${selected_pheno}_${TEST_POP}_PGS.txt"
    pgs_valid_thresh="${best_thresh}${selected_pheno}_thresh.txt"
    best_pgs_thresh=$(head -2 "$processed_pgs" | tail -1 | awk '{print $2}' | tr -d '"')
    echo "$best_pgs_thresh 0 $best_pgs_thresh" > "$pgs_valid_thresh"
    echo "$best_pgs_thresh"  

    # build pgs
    pgs_out_file="${pgs_path}${selected_pheno}_${valid_pop}_pgs_valid"
    last_thresh="${THRESH[-1]}" #grab ex. thresh
    gwas_file="${JOINED_GWAS}${selected_pheno}/chr1-22_${TRAIN_POP}_${selected_pheno}.assoc.linear"
    valid_snp_file="${CLUMPED_GENOS}ukb_chr1-22_${TEST_POP}_${selected_pheno}.valid.snp" #use test pgs clumped snps
    
    # for the --score flag : 3 5 12 grabs snp ID, alt allele, & effect size estimate
    # for --q-score-range : cols 3 & 15 of "$gwas_file" grabs snp IDs & pval cols
    if [[ ! -f "$pgs_out_file"."$best_pgs_thresh".sscore ]]; then
      ${LOCAL_PLINK2} \
        --bfile ${pop_path}${pop_file}\
        --keep ${genos}${pop_file}.fam \
        --score "$gwas_file" 3 5 12 header cols=+scoresums \
        --q-score-range "$pgs_valid_thresh" "$gwas_file" 3 15 header \
        --extract "$valid_snp_file" \
        --out "$pgs_out_file"
      echo "Std PGS complete"
      #add back to cut pop to 1K
      #--keep ${pop_path}/${pop_file}_QC_1000.fam \
    else
      echo "std pgs already exists"
    fi
    
  done
  echo "valid PGS done"



# Calculate total runtime in seconds
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds" 