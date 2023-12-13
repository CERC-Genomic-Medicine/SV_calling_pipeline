#!/usr/bin/env nextflow 

process FiltTest {
   
   cache "lenient"
   cpus 1
   
   
   publishDir "${params.result}/prepareVCFsToMerge"

   input:   
   val(input_label)

   output:
   tuple path("${input_label}.*.vcf*"), val(input_label)

   script:
   """
      tools=("delly" "manta")
      variants=("BND" "DEL" "DUP" "INS" "INV")

      for tool in \${tools[@]};
      do
         bcftools view ${params.result}/\${tool}/${input_label}.\${tool}.vcf.gz -r chr1, ... chr22, chrX, chrY > ${input_label}.\${tool}.vcf.01 
         ${params.survivor} filter ${input_label}.\${tool}.vcf.01 NA 50 -1 0 -1 ${input_label}.\${tool}.vcf.02
         awk -F '\t' '{if(\$0 ~ "#") print; else if(\$7 == "PASS") print}' ${input_label}.\${tool}.vcf.02 > ${input_label}.\${tool}.vcf.03
         Rscript ${params.scripts}/changeSampleName.R ${input_label}.\${tool}.vcf.03 \${tool}_${input_label} ${input_label}.\${tool}.vcf.04
         Rscript ${params.scripts}/renameRefAlt.R ${input_label}.\${tool}.vcf.04 ${input_label}.\${tool}.vcf.05.gz
         gunzip ${input_label}.\${tool}.vcf.05.gz
         cp ${input_label}.\${tool}.vcf.05 ${input_label}.\${tool}.vcf.filt
         bgzip -c ${input_label}.\${tool}.vcf.filt > ${input_label}.\${tool}.vcf.filt.gz
         tabix -p vcf ${input_label}.\${tool}.vcf.filt.gz
         bcftools query -f '%CHROM %POS %INFO/SVTYPE [%GT]\n' ${input_label}.\${tool}.vcf.filt.gz > ${input_label}.\${tool}.vcf.filt.genotypes
         sed 's/SVTYPE=DUP/SVTYPE=INS/' ${input_label}.\${tool}.vcf.filt | grep -E '^#|SVTYPE=INS' > ${input_label}.\${tool}.vcf.filt.DUPtoINS
            for val in \${variants[@]}; 
            do
              grep -E '^#|SVTYPE=\${val}' ${input_label}.\${tool}.vcf.filt > ${input_label}.\${tool}.vcf.filt.\${val}
            done
       done
   """
}
  
process MergeSVs {

   input:   
   val(input_label)

    output: 
    path("*.callsMerged_*.vcf")

    publishDir "${params.result}/prepareVCFsToMerge/SURVIVOR"

    script:
    """
    tools=("delly" "manta")
    variants=("BND" "DEL" "DUP" "INS" "INV")

    for tool in \${tools[@]};
    do
      find ${params.result}/prepareVCFsToMerge/ -name "${input_label}.\${tool}.vcf.filt*" -type f > vcfCallFiles_ALL.list
      ${params.survivor} merge vcfCallFiles_ALL.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} callsMerged_ALL.vcf
      vcf-sort -c callsMerged_ALL.vcf > callsMerged_ALL.sorted.vcf
      bgzip -c callsMerged_ALL.sorted.vcf > \${tool}.callsMerged_ALL.sorted.vcf.gz
      tabix -p vcf \${tool}.callsMerged_ALL.sorted.vcf.gz
      for var in \${variants[@]} 
      do
         find ${params.result}/prepareVCFsToMerge/ -name "${input_label}.*.vcf.filt.\${var}" -type f > vcfCallFiles_\${var}.list
         ${params.survivor} merge vcfCallFiles_\${var}.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} \${tool}.callsMerged_\${var}.vcf
      done
    done
    """
    }


   //  bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO callsMerged_ALL.sorted.vcf.gz > callsMerged_ALL.sorted.clean.vcf

   //  bgzip -c callsMerged_ALL.sorted.clean.vcf > \${tool}.callsMerged_ALL.sorted.clean.vcf.gz
   //  tabix -p vcf \${tool}.callsMerged_ALL.sorted.clean.vcf.gz


   workflow {

   // inputFiles = Channel.fromPath(params.inputFiles).map{file -> file.getSimpleName()}
   inputFiles = Channel.fromPath(params.inputFilesSubset).map{file -> file.getSimpleName()}
   inputFiles | FiltTest | MergeSVs
}
