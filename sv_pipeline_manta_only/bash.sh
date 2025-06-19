#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=12:00:00
#SBATCH --account=rrg-vmooser
module load nextflow
module load bcftools
module load java
module load r
module load python
module load tabix 

NXF_DISABLE_CHECK_LATEST=true nextflow run manta_filtering_post_gt.nf -resume 