#!/usr/bin/env nextflow 

// Before running this script load samtools:
// module load samtools
// module load bcftools
// module load delly
// module load gcc
// module load manta
// line 62 # can you stream output to stdout


// process CramToBam {
//    errorStrategy 'retry'
//    maxRetries 3
//    cache "lenient"
//    cpus 1
//    memory "16 GB"
//    time "24h"
   
//    input:   
//    tuple val(input_label), path(cram), path(crai)

//    output:
//    tuple val(input_label), path("*.bam"), path("*.bam.bai") 
   
//    script:
//    """
//    samtools view -b -T ${params.referenceGenome} ${cram} > ${params.bamsFolder}/${cram.getBaseName()}.bam
//    samtools index -b ${cram.getBaseName()}.bam
//    """
// }

process DellySV {
   maxRetries 0
   cache "lenient"
   cpus 16
   memory "64 GB"
   time "48h"

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   path("${input_label}.delly.vcf.gz*")

   publishDir "${params.result}/delly", pattern: '*.vcf.gz*', mode: 'copy'

   script:    
   """

   export OMP_NUM_THREADS=16

   ${params.delly} call -o ${input_label}.bcf -g ${params.referenceGenome} ${cram}  
   bcftools view ${input_label}.bcf -Oz -o ${input_label}.delly.vcf.gz
   bcftools index -t ${input_label}.delly.vcf.gz
   """
} 

process MantaSV {

   cache "lenient"
   cpus 16
   memory "64 GB"
   time "48h"


   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   path("${input_label}.manta.vcf.gz*")

   publishDir "${params.result}/manta" 

   """

   configManta.py --bam ${cram} --referenceFasta ${params.referenceGenome}
   ./MantaWorkflow/runWorkflow.py -j 4
   mv ./MantaWorkflow/results/variants/diploidSV.vcf.gz ${input_label}.manta.vcf.gz
   
   """
}
   // bcftools index -t ${input_label}.manta.vcf.gz
workflow {
   
   inputFiles = Channel.fromPath(params.inputFiles).map{file -> [file.getSimpleName(), file, file + ".crai"]}
// inputFiles.view()

   inputFiles | MantaSV
   inputFiles | DellySV
}
