#!/usr/bin/env nextflow 

// Before running this script load samtools:
// module load samtools
// module load bcftools
// module load delly
// module load gcc
// module load manta


process DellySV {
// call Delly, compress and index resulting file 
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "8h"

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val("${input_label}"), file("${input_label}.delly.vcf.gz"), file("${input_label}.delly.vcf.gz.tbi")

   publishDir "${params.result}/dellySampleSubset", pattern: '*.vcf.gz*', mode: 'copy'

   script:    
   """
   ${params.delly} call -o ${input_label}.bcf -g ${params.referenceGenome} ${cram}  
   bcftools view ${input_label}.bcf -Oz -o ${input_label}.delly.vcf.gz
   bcftools index -t ${input_label}.delly.vcf.gz
   """
} 

process FilterDellySV {
   cache "lenient"
   cpus 1
   memory "4 GB"
   time "2h"

   input:
   tuple val(input_label), file(delly_sv_vcf), file(delly_sv_vcf_tbi)
   
   publishDir "${params.result}/FilterSampleSubset"

   output:
   tuple path("${input_label}.delly.vcf.filt*"), val(input_label)

   //filter the vcf file to remove any SVs smaller than 50 bp and any SVs without 'PASS'
   // separate the vcf file into different files by SVTYPE

   script:
   """
      variants=("BND" "DEL" "DUP" "INS" "INV")

        bcftools view ${delly_sv_vcf} -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${input_label}.delly01.vcf      
        ${params.survivor} filter ${input_label}.delly01.vcf NA 50 -1 0 -1 ${input_label}.delly02.vcf
        awk -F '\t' '{if(\$0 ~ "#") print; else if(\$7 == "PASS") print}' ${input_label}.delly02.vcf > ${input_label}.delly03.vcf
        bcftools query -l ${input_label}.delly03.vcf | awk '{ print \$1, "'"delly"'_"\$1 }' > new_sample_names.txt
        bcftools reheader --samples new_sample_names.txt ${input_label}.delly03.vcf > ${input_label}.delly04.vcf


        sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.delly04.vcf | grep -E '^#|SVTYPE=INS' > ${input_label}.delly.vcf.filt.DUPtoINS
         
        for val in \${variants[@]}; 
        do
            grep -E '^#|SVTYPE='"\${val}" ${input_label}.delly04.vcf > ${input_label}.delly.vcf.filt.\${val}
        done
   """

}


process MantaSV {

   cache "lenient"
   cpus 4
   memory "4 GB"
   time "8h"


   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.manta.vcf.gz"), path("${input_label}.manta.vcf.gz.tbi")

   publishDir "${params.result}/mantaSampleSubset" 

   """

   configManta.py --bam ${cram} --referenceFasta ${params.referenceGenome}
   ./MantaWorkflow/runWorkflow.py -j 4
   mv ./MantaWorkflow/results/variants/diploidSV.vcf.gz ${input_label}.manta.vcf.gz
   bcftools index -t ${input_label}.manta.vcf.gz 
   """
}

process FilterMantaSV {
// filter the vcf file to remove any SVs smaller than 50 bp and any SVs without 'PASS'
// separate the vcf file into different files by SVTYPE

   cache "lenient"
   cpus 1
   memory "4 GB"
   time "2h"

   input:
   tuple val(input_label), path(manta_sv_vcf), path(manta_sv_vcf_tbi)
   
   publishDir "${params.result}/FilterSampleSubset"

   output:
   tuple path("${input_label}.manta.vcf.filt*"), val(input_label)

   script:
   """
      variants=("BND" "DEL" "DUP" "INS" "INV")

        bcftools view ${manta_sv_vcf} -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${input_label}.manta01.vcf      
        ${params.survivor} filter ${input_label}.manta01.vcf NA 50 -1 0 -1 ${input_label}.manta02.vcf
        awk -F '\t' '{if(\$0 ~ "#") print; else if(\$7 == "PASS") print}' ${input_label}.manta02.vcf > ${input_label}.manta03.vcf
        bcftools query -l ${input_label}.manta03.vcf | awk '{ print \$1, "'"manta"'_"\$1 }' > new_sample_namesmanta.txt
        bcftools reheader --samples new_sample_namesmanta.txt ${input_label}.manta03.vcf > ${input_label}.manta04.vcf


        sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.manta04.vcf | grep -E '^#|SVTYPE=INS' > ${input_label}.manta.vcf.filt.DUPtoINS
         
        for val in \${variants[@]}; 
        do
            grep -E '^#|SVTYPE='"\${val}" ${input_label}.manta04.vcf > ${input_label}.manta.vcf.filt.\${val}
        done
   """
}



workflow {
   
   inputFiles = Channel.fromPath(params.inputFilesSubset).map{file -> [file.getSimpleName(), file, file + ".crai"]}

   inputFiles | MantaSV | FilterMantaSV 
   inputFiles | DellySV | FilterDellySV 
}

