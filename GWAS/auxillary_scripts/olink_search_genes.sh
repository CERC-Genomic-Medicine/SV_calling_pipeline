#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=72:00:00
#SBATCH --account=rrg-vmooser

# last updated: Feb 28 2025 


module load bcftools
olink_csv=/home/oliviagc/projects/rrg-vmooser/CERC_Private/Omics/BQC19/OLINK/Olink_T48/BQC19_2_NPX.csv

while IFS=$'\t' read gene; do
    echo ${gene}
    awk -F';' '$5 ~ /'${gene}'/ {print}' ${olink_csv}
done < /home/oliviagc/scratch/genes.txt
