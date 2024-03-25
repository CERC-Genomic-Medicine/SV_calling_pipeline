#!/usr/bin/env nextflow 

process DellySV {
// call Delly, compress and index resulting file 
   cache "lenient"
   cpus 6
   memory "20 GB"
   time "48h"

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val("${input_label}"), file("${input_label}.delly.vcf.gz"), file("${input_label}.delly.vcf.gz.tbi")

   publishDir "${params.result}/delly", pattern: '*.vcf.gz*', mode: 'copy'

   script:    
   """
   ${params.delly} call -o ${input_label}.bcf -g ${params.referenceGenome} ${cram}  
   bcftools view ${input_label}.bcf -Oz -o ${input_label}.delly.vcf.gz
   bcftools index -t ${input_label}.delly.vcf.gz
   """
} 


process MantaSV {

   cache "lenient"
   cpus 4
   memory "8 GB"
   time "12h"


   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.manta.vcf.gz"), path("${input_label}.manta.vcf.gz.tbi")

   publishDir "${params.result}/manta", pattern: '*.vcf.gz*', mode: 'copy'

   """
   # call manta, restricted to ref regions included in bed file 
   ${params.manta}/bin/configManta.py --bam ${cram} --referenceFasta ${params.referenceGenome} --callRegions ${params.referenceDir}/GRCh38.include.bed.gz
   # execute script with specified job count equal to 4 
   ./runWorkflow.py -j 4
   
   # analysis of BNDs to determine INVs, output = vcf
   ${params.manta}/libexec/convertInversion.py ${params.manta}/libexec/samtools ${params.referenceGenome} ./MantaWorkflow/results/variants/diploidSV.vcf.gz > ${input_label}.manta.vcf
   
   #compression to vcf.gz 
   bcftools view ${input_label}.manta.vcf -Oz -o ${input_label}.manta.vcf.gz
   
   #tabix indexation
   bcftools index -t ${input_label}.manta.vcf.gz
   """
}




 process LumpySV {

   cache "lenient"
   cpus 4
   memory "16 GB"
   time "12h"
   scratch '$SLURM_TMPDIR'
  

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.lumpy.vcf.gz"), path("${input_label}.lumpy.vcf.gz.tbi")

   publishDir "${params.result}/lumpy", pattern: '*.vcf.gz*', mode: 'copy'

   script: 
   """
   chrom=chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}

   samtools view -h -b -F 1294 ${cram} \${chrom} > ${input_label}.discordant.unsorted.bam
   
   samtools view -h ${cram} \${chrom} | ${params.lumpy}/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${input_label}.splitters.unsorted.bam
   
   samtools sort -o ${input_label}.discordant.bam ${input_label}.discordant.unsorted.bam
   samtools sort -o ${input_label}.splitters.bam ${input_label}.splitters.unsorted.bam 

   sed "s|SAMTOOLS=|SAMTOOLS=\$(which samtools)|" ${params.lumpy}/bin/lumpyexpress.config > lumpyexpress.config
   sed -i "s|PYTHON=.*|PYTHON=\$(which python)|" lumpyexpress.config
   
   ${params.lumpy}/bin/lumpyexpress -B ${cram} -S ${input_label}.discordant.bam -D ${input_label}.splitters.bam -o ${input_label}.lumpy.vcf -K lumpyexpress.config 
   bgzip ${input_label}.lumpy.vcf
   tabix -p vcf ${input_label}.lumpy.vcf.gz
   """
 }


workflow {
   
   inputFiles = Channel.fromPath(params.inputFilesSubsetFourteen).map{file -> [file.getSimpleName(), file, file + ".crai"]}
   inputFiles | MantaSV 
   inputFiles | DellySV 
   inputFiles | LumpySV 
}