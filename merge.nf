#!/usr/bin/env nextflow
 process MergeTools {
   // first we take all the SV types for each sample and put them in a list by SV type 
  //  cache "lenient"
  //  cpus 2
  //  memory "8 GB"
  //  time "8h"

    input:   
    tuple val(input_label), val(variant), path(vcf_file)

    output:
    tuple val(variant), path("${input_label}.${variant}.merge.fixed.vcf")
   
   publishDir "${params.result}/mergeToolsTest"
     script:
     """
      ls ${vcf_file} > ${input_label}.${variant}.vcf.list
      ${params.survivor} merge ${input_label}.${variant}.vcf.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} ${input_label}.${variant}.merge.vcf

      bcftools sort -Oz -o ${input_label}.${variant}.merge.sorted.vcf.gz ${input_label}.${variant}.merge.vcf 
      bcftools index -t ${input_label}.${variant}.merge.sorted.vcf.gz
      Rscript ${params.scripts}/fixSURVIVOR.R ${input_label}.${variant}.merge.sorted.vcf.gz ${input_label} ${input_label}.${variant}.merge.fixed.vcf.gz
      gunzip -f ${input_label}.${variant}.merge.fixed.vcf.gz
    """
}


process MergeSamples {
    // cache "lenient"
    // cpus 1
    // memory "4 GB"
    // time "8 h"
    
    input:   
    tuple val(variant), path(merged_tools)

    output:
    path("${variant}.samples_merged.vcf")

    publishDir "${params.result}/mergeSamples"
    script:
    """
      ls ${merged_tools} > ${variant}.allsamples.vcf.list
      ${params.survivor} merge ${variant}.allsamples.vcf.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} 0 ${variant}.samples_merged.vcf 
    """
}


workflow {
inputFiles = Channel
                .fromPath(params.inputFilesMerge)
                .map{file -> [file.getSimpleName(), file.getExtension(), file]}
                .groupTuple(by:[0,1] )
                
inputFiles | MergeTools | groupTuple(by:0) | MergeSamples 
}