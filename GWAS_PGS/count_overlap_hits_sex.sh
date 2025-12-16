#!/usr/bin/env bash

#SBATCH -J compare_gwas
#SBATCH --time=20:00:00
#SBATCH --mem-per-cpu=75GB
#SBATCH --partition=tier2q
#SBATCH -o /home/mcostantino/job-output/run_PGS/compare_gwas_loop.out
#SBATCH -e /home/mcostantino/job-output/run_PGS/compare_gwas_loop.err

pheno="testosterone"

LOCAL_PLINK2="/ess/home/home1/mcostantino/plink2"
threshold=5e-8
output_file="/gpfs/data/ukb-share/dahl/manuela/scale_project/overlap_testo_female_log.txt"

# Clear previous output
> "$output_file"


echo "Processing phenotype: $pheno"

gwas1="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_female/${pheno}/chr1-22_whitebrit_${pheno}_females.glm.linear"
gwas2="/gpfs/data/ukb-share/dahl/manuela/joined_gwas_results_log/${pheno}/chr1-22_whitebrit_${pheno}.assoc.linear"
union_gwas="/scratch/mcostantino/descaler/union_gwas/${pheno}_female_default.assoc.linear"
clump_file="/scratch/mcostantino/descaler/union_gwas/clumped_genos/${pheno}_female_default"
valid_snp_file="/scratch/mcostantino/descaler/union_gwas/clumped_genos/${pheno}_female_default.valid.snp"

mkdir -p "$(dirname "$union_gwas")"
mkdir -p "$(dirname "$clump_file")"

# Take min pvalue of each
awk 'BEGIN{OFS="\t"} NR==FNR{a[$3]=$15;next} ($3 in a){$15=($15<a[$3]?$15:a[$3]); print}' "$gwas1" "$gwas2" > "$union_gwas"

# Clump file if it doesn't exist
if [[ ! -f "$clump_file".clumps ]]; then
    ${LOCAL_PLINK2} \
        --bfile /gpfs/data/ukb-share/genotypes/pop_genos/whitebrit/ukb_chr1-22_whitebrit_QC \
        --clump-p1 1 \
        --clump-r2 0.1 \
        --clump-kb 250 \
        --clump "$union_gwas" \
        --out "$clump_file"
    echo "$pheno clump done"
    awk 'NR!=1{print $3}' ${clump_file}.clumps > "$valid_snp_file"
else
    echo "Output file already exists at $clump_file"
fi

# Subset to clumped SNPs
awk 'NR==FNR {keep[$1]; next} ($2 in keep)' "$valid_snp_file" "$union_gwas" > "${union_gwas%.assoc.linear}_clumped.assoc.linear"
awk 'NR==FNR {keep[$1]; next} ($2 in keep)' "$valid_snp_file" "$gwas1" > "${gwas1%.assoc.linear}_clumped.assoc.linear"
awk 'NR==FNR {keep[$1]; next} ($2 in keep)' "$valid_snp_file" "$gwas2" > "${gwas2%.assoc.linear}_clumped.assoc.linear"

# Identify significant SNPs
awk -v t=$threshold 'NR==FNR {snp[$1]; next} ($3 in snp) && ($15 < t) {print $3}' "$valid_snp_file" "$union_gwas" > sig_union.tmp
awk -v t=$threshold 'NR==FNR {snp[$1]; next} ($3 in snp) && ($15 < t) {print $3}' "$valid_snp_file" "$gwas1" > sig_gwas1.tmp
awk -v t=$threshold 'NR==FNR {snp[$1]; next} ($3 in snp) && ($15 < t) {print $3}' "$valid_snp_file" "$gwas2" > sig_gwas2.tmp

sig_union=$(wc -l < sig_union.tmp)
sig_gwas1=$(wc -l < sig_gwas1.tmp)
sig_gwas2=$(wc -l < sig_gwas2.tmp)

# Overlap counts
shared_1_union=$(grep -Fxf sig_gwas1.tmp sig_union.tmp | wc -l)
shared_2_union=$(grep -Fxf sig_gwas2.tmp sig_union.tmp | wc -l)
shared_1_2=$(grep -Fxf sig_gwas1.tmp sig_gwas2.tmp | wc -l)
shared_all=$(grep -Fxf sig_gwas1.tmp sig_gwas2.tmp | grep -Fxf - sig_union.tmp | wc -l)

# Append to output
echo -e "${pheno}\t${shared_1_union}\t${shared_2_union}\t${shared_all}" >> "$output_file"

# Cleanup
rm -f sig_union.tmp sig_gwas1.tmp sig_gwas2.tmp
