#!/usr/bin/env bash

#SBATCH -J std_gwas
#SBATCH --mem=10GB
#SBATCH --time=08:00:00
#SBATCH --partition=tier1q
#SBATCH --array=1-22%3
#SBATCH -o /home/mcostantino/gwas/std_gwas_%a_%A.out
#SBATCH -e /home/mcostantino/gwas/std_gwas_%a_%A.err

now=$(date)
# Record job start time
start_time=$(date+%s)

echo "Date run: $now"

lambda=0
echo "Lambda: $lambda"

pheno=$1
GENOS="/scratch/reneefonseca/genotypes/"
GWAS="/gpfs/data/ukb-share/dahl/manuela/gwas_results/"
GWAS_LOG="/gpfs/data/ukb-share/dahl/manuela/gwas_results_log/"
LOCAL_PLINK2="/ess/home/home1/mcostantino/plink2"
TRAIN_POP="whitebrit"
INTERMEDIATE_DIR="/scratch/mcostantino/descaler/intermediate_files"
chr=${SLURM_ARRAY_TASK_ID:-21}
pop_path="${GENOS}${TRAIN_POP}/"
pop_file="ukb_chr${chr}_${TRAIN_POP}_QC"
keep_file="/gpfs/data/ukb-share/genotypes/pop_genos/${TRAIN_POP}/ukb_chr1-22_${TRAIN_POP}_QC_train"

echo "Processing pheno ${pheno} and chr ${chr}"

# Function to process phenotype name
PROCESS_PHENO_NAME() {
    local pheno_name=$1
    local phenoNoDigits=${pheno_name%%[0-9]*}
    local phenoLower=$(echo "$phenoNoDigits" | tr '[:upper:]' '[:lower:]')

    case $pheno_name in
        *"FEV1674178"*) phenoNoDigits="FEV1"; phenoLower="fev1" ;;
        *"IGF-1674178"*) phenoNoDigits="IGF-1"; phenoLower="igf-1" ;;
        *"HbA1c674178"*) phenoNoDigits="HbA1c"; phenoLower="hba1c" ;;
        *"EA4"*) phenoNoDigits="EA4"; phenoLower="ea4" ;;
        *"A1c674178"*) phenoNoDigits="A1c"; phenoLower="a1c" ;;
    esac

    echo "$phenoNoDigits $phenoLower"
}

# Process the selected phenotype's name
parse_pheno=$(PROCESS_PHENO_NAME "$pheno")
phenoNoDigits=${parse_pheno%% *}
phenoLower=${parse_pheno##* }

# Get phenotype file
pheno_file="/gpfs/data/ukb-share/extracted_phenotypes/${phenoNoDigits}/${pheno}.pheno"
echo "Phenotype file: $pheno_file"

# Get phenotype code
pheno_code=$(head -1 "$pheno_file" | awk '{ print $3 }')

# Apply log scale if lambda is 0
if (( $(echo "$lambda == 0" | bc -l) )); then
    echo "Lambda is 0. Generating log-scaled phenotype file."
    mkdir -p "$INTERMEDIATE_DIR"
    pheno_file_log="${INTERMEDIATE_DIR}/${pheno}_log.pheno"

    # Log-transform phenotype values
    awk '
    BEGIN { OFS="\t" }
    NR == 1 { print $0 }  # Keep header row
    NR > 1 {
        if ($3 ~ /^[0-9.]+$/ && $3 > 0) {
            $3 = log($3)
        } else {
            $3 = "NA"  # Replace invalid entries
        }
        print $0
    }' "$pheno_file" > "$pheno_file_log"

    pheno_file="$pheno_file_log"
    echo "Log-scaled phenotype file saved to: $pheno_file"
fi


# Define GWAS output path based on lambda
if (( $(echo "$lambda == 0" | bc -l) )); then
    gwas_path="${GWAS_LOG}${phenoLower}/"
else
    gwas_path="${GWAS}${phenoLower}/"
fi
mkdir -p "$gwas_path"

# Covariates file
covar_file="/gpfs/data/ukb-share/extracted_phenotypes/covariates_sa40PC/covariates_sa40PC674178.pheno"
ncovar='3-14'

# Define the output file path
out_file="${gwas_path}chr${chr}_${TRAIN_POP}_${phenoLower}"

# Run PLINK2
if [[ ! -f "$out_file.$pheno_code.glm.linear" ]]; then
    ${LOCAL_PLINK2} --bfile ${pop_path}${pop_file} \
        --keep ${keep_file} \
        --pheno "$pheno_file" \
        --pheno-name "$pheno_code" \
        --no-input-missing-phenotype \
        --glm hide-covar \
        --covar-variance-standardize 34-0.0 \
        --covar "$covar_file" \
        --covar-col-nums "$ncovar" \
        --out "$out_file"
else
    echo "Output file already exists: $out_file.$pheno_code.glm.linear"
fi

# Calculate total runtime
end_time=$(date +%s)
runtime=$((end_time - start_time))
echo "Job Runtime: $runtime seconds"