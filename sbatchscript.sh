#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G
#SBATCH --time=120:00:00
#SBATCH --account=def-vmooser

module load StdEnv/2020
module load gcc/9.3.0
module load perl/5.30.2
module load bcftools
module load samtools
module load root/6.26.06
module load bowtie2
module load java
module load nextflow
module load apptainer

source /home/oliviagc/projects/rrg-vmooser/oliviagc/venv_py2.7.18/bin/activate

nextflow run pipeline.nf -with-report pipeline_alltools_batch2.html
