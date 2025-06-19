#!/bin/bash

module load bcftools
for file in /home/oliviagc/scratch/samples_merged_*.Final.remove_cpx.vcf.gz #change path to your vcf input files 
do
    filename=$(basename "$file" .vcf.gz)                     
    bcftools annotate --rename-chrs chr_name_conv.txt "$file" -Oz -o "/home/oliviagc/scratch/BQC/BQC19_SAIGE/inputs_saige_all/input_vcfs_SVs/${filename}.vcf.gz"  #change to the folder you want to put your new vcfs in 
    bcftools index "/home/oliviagc/scratch/BQC/BQC19_SAIGE/inputs_saige_all/input_vcfs_SVs/${filename}.vcf.gz"  #change to the folder you want to put your new vcfs in 
done