#!/bin/bash

module load bcftools
for file in /home/oliviagc/scratch/BQC_SAIGE/BQC_SAIGE_vcfs/*Final.gd_overlap.vcf.gz #change path to your vcf input files 
do
    filename=$(basename "$file" .vcf.gz)   
    for chr in {1..22} X    
    do               
        bcftools view -i "CHROM = '$chr'" "$file" -Oz -o "./input_vcfs_oc/${filename}.${chr}.vcf.gz"  #change to the folder you want to put your new vcfs in 
        bcftools index "./input_vcfs_oc/${filename}.${chr}.vcf.gz"  #change to the folder you want to put your new vcfs in 
    done                  
done