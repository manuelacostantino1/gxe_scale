#!/usr/bin/env bash

#SBATCH -J join_gwas
#SBATCH --mem=20GB
#SBATCH --time=10:00:00
#SBATCH --partition=tier1q
#SBATCH -o /home/mcostantino/join_gwas/join_gwas_%A.out
#SBATCH -e /home/mcostantino/join_gwas/join_gwas_%A.err

now=$(date)
# Record job start time
start_time=$(date +%s)

echo "Date run: $now"

lambda=1
echo "Lambda: $lambda"

pheno=$1
GWAS="/gpfs/data/ukb-share/dahl/manuela/gwas_results/"
GWAS_LOG="/gpfs/data/ukb-share/dahl/manuela/gwas_results_log/"
LOCAL_PLINK2="/ess/home/home1/mcostantino/plink2"
TRAIN_POP="whitebrit"
JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results/"

if (( $(echo "$lambda == 0" | bc -l) )); then
    gwas_path="${GWAS_LOG}${phenoLower}/"
	GWAS="$GWAS_LOG"
	JOINED_GWAS="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results_log/"
else
    gwas_path="${GWAS}${phenoLower}/"
fi

join_gwas(){	
	# set up directories
	input_dir="${GWAS}${pheno}/"
	output_dir="${JOINED_GWAS}${pheno}/"
	mkdir -p "$output_dir"

	output_file="${output_dir}chr1-22_${TRAIN_POP}_${pheno}.assoc.linear"
	echo "$output_file"
	if [ ! -f "$output_file" ]; then
		first_file=$(ls ${input_dir}*${context}*.linear | head -1) 
		echo "$first_file"
		file_end=$(echo "$first_file" | cut -d'.' -f2-)
		chr1_gwas="${input_dir}chr1_${TRAIN_POP}_${pheno}.${file_end}"
		
		echo "joining ${chr1_gwas} and all other chromosomes"
		# write header to output file
		head -1 "$chr1_gwas" > "$output_file"
		# append chr data to output file
		tail -n +2 -q "${input_dir}"chr{1..22}_"${TRAIN_POP}"_"${pheno}"."${file_end}" >> "$output_file"
	else
		echo "${pheno} join gwas already made"
	fi
}

join_gwas "$pheno"