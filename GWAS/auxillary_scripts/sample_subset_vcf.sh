#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --account=rrg-vmooser
#SBATCH --mem-per-cpu=8G
set -e

module load bcftools


sample_subset=/home/oliviagc/scratch/SV_calling_pipeline/BQC_SAIGE/BQC19_SAIGE/sample_list_oc/bqc_f

for sv in DEL BND INS INV DUP
do
    VCF=/home/oliviagc/scratch/SV_calling_pipeline/BQC_SAIGE/BQC19_SAIGE/input_vcfs/samples_merged_${sv}.Final.gd_overlap.vcf.gz
    OUT=/home/oliviagc/scratch/SV_calling_pipeline/BQC_SAIGE/BQC19_SAIGE/input_vcfs_female/samples_merged_${sv}.Final.gd_overlap_female.vcf.gz
    bcftools view -S ${sample_subset} ${VCF} -Oz -o ${OUT}
    bcftools index ${OUT}
done
