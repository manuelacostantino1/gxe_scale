#!/usr/bin/env bash

#SBATCH -J t0
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --partition=tier1q
#SBATCH -o /home/mcostantino/join_gwas/join_gwas_%A.out
#SBATCH -e /home/mcostantino/join_gwas/join_gwas_%A.err

now=$(date)
start_time=$(date +%s)
echo "Date run: $now"

pheno=$1
TRAIN_POP="whitebrit"

# Input GWAS dirs
GWAS_MALE="/gpfs/data/ukb-share/dahl/manuela/gwas_results_male/"
GWAS_FEMALE="/gpfs/data/ukb-share/dahl/manuela/gwas_results_female/"
GWAS_LOG_MALE="/gpfs/data/ukb-share/dahl/manuela/gwas_results_log_male/"
GWAS_LOG_FEMALE="/gpfs/data/ukb-share/dahl/manuela/gwas_results_log_female/"

# Output joined dirs
JOINED_GWAS_MALE="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_male/"
JOINED_GWAS_FEMALE="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_female/"
JOINED_GWAS_MALE_LOG="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_male_log/"
JOINED_GWAS_FEMALE_LOG="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_female_log/"

# Methods: default and log
for method in default log; do
  for sex in male female; do
    
    # Select input/output directories
    if [[ "$method" == "default" && "$sex" == "male" ]]; then
      input_dir="${GWAS_MALE}${pheno}/"
      output_dir="${JOINED_GWAS_MALE}${pheno}/"
      suffix="males"   # adjust if your default files actually lack "_males"
    elif [[ "$method" == "default" && "$sex" == "female" ]]; then
      input_dir="${GWAS_FEMALE}${pheno}/"
      output_dir="${JOINED_GWAS_FEMALE}${pheno}/"
      suffix="females"
    elif [[ "$method" == "log" && "$sex" == "male" ]]; then
      input_dir="${GWAS_LOG_MALE}${pheno}/"
      output_dir="${JOINED_GWAS_MALE_LOG}${pheno}/"
      suffix="males"
    else
      input_dir="${GWAS_LOG_FEMALE}${pheno}/"
      output_dir="${JOINED_GWAS_FEMALE_LOG}${pheno}/"
      suffix="females"
    fi

    mkdir -p "$output_dir"
    output_file="${output_dir}chr1-22_${TRAIN_POP}_${pheno}_${suffix}.glm.linear"

    if [[ ! -f "$output_file" ]]; then
      # detect file ending (everything after <pheno>_<suffix>.)
      first_file=$(ls ${input_dir}chr1_${TRAIN_POP}_${pheno}_${suffix}*.glm.linear 2>/dev/null | head -1)

      if [[ -z "$first_file" ]]; then
        echo "⚠ No files found in ${input_dir} for ${sex}, ${method}!"
        continue
      fi

      file_end=$(echo "$first_file" | sed -E "s/.*${pheno}_${suffix}\.//")

      chr1_gwas="${input_dir}chr1_${TRAIN_POP}_${pheno}_${suffix}.${file_end}"
      echo "Joining ${chr1_gwas} and other chromosomes for ${sex} (${method})..."

      # write header
      head -1 "$chr1_gwas" > "$output_file"
      # append all chromosomes
      for chr in {1..22}; do
        tail -n +2 "${input_dir}chr${chr}_${TRAIN_POP}_${pheno}_${suffix}.${file_end}" >> "$output_file"
      done
    else
      echo "${pheno} ${sex} ${method} joined GWAS already exists"
    fi

  done
done

end_time=$(date +%s)
echo "Elapsed time: $(( (end_time - start_time)/60 )) minutes"