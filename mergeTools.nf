#!/usr/bin/env nextflow

process MergeTool {
   cache "lenient"
   cpus 1
   input:   
   val(input_label)

   output:
   path("${input_label}.*.merged.sorted.vcf3.gz")
   
   publishDir "${params.result}/mergeTools"

   script:
   """
    variants=("BND" "DEL" "DUP" "INS" "INV")

      for var in \${variants[@]}; 
      do
         find ${params.result}/prepareVCFsToMerge -type f -name "${input_label}.*.vcf.filt.\${var}" > ${input_label}.filtered.vcf.list1.\${var}
         ${params.survivor} merge ${input_label}.filtered.vcf.list1.\${var} ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} ${input_label}.merged.vcf1.\${var}
         sed -i 's|FORMAT=<ID=DR,Number=1,Type=Integer|FORMAT=<ID=DR,Number=1,Type=String|g' ${input_label}.merged.vcf1.\${var}
         sed -i 's|ID=LN,Number=1,Type=Integer|ID=LN,Number=1,Type=String|g' ${input_label}.merged.vcf1.\${var}
         
         find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.*.vcf.filt.\${var}" > ${input_label}.filtered.vcf.list2.\${var}
         ${params.survivor} merge ${input_label}.filtered.vcf.list2.\${var} ${params.breakpoint_dist} 2 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} ${input_label}.\${var}.merged.vcf2
        
         find . -type f -name "${input_label}.filtered.vcf.list[1|2].\${var}" > ${input_label}.filtered.vcf.list.\${var}
         ${params.survivor} merge ${input_label}.filtered.vcf.list.\${var} ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} ${input_label}.\${var}.merged.vcf3
         vcf-sort -c ${input_label}.\${var}.merged.vcf3 > ${input_label}.\${var}.merged.sorted.vcf3
         bgzip -f -c ${input_label}.\${var}.merged.sorted.vcf3 > ${input_label}.\${var}.merged.sorted.vcf3.gz
         tabix -p vcf ${input_label}.\${var}.merged.sorted.vcf3.gz
      done 
   """     
}

      //    Rscript ${params.scripts}/fixSURVIVOR.R ${input_label}.\${var}.merged.sorted.vcf3.gz ${input_label} ${input_label}.\${var}.merged.vcf.gz
      //    gunzip -f ${input_label}.\${var}.merged.vcf.gz

workflow {
// inputFiles = Channel.fromPath(params.inputFiles).map{file -> file.getSimpleName()}
inputFiles = Channel.fromPath(params.inputFilesSubset).map{file -> file.getSimpleName()}
inputFiles | MergeTool
 
}