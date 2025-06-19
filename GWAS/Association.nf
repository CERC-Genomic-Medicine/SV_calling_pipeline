#!/usr/bin/env nextflow

// module load nextflow apptainer 

// formatic snp array file into plink/bed format
process PlinkFormatting {
  cache "lenient"
  cpus 1
  scratch '$SLURM_TMPDIR'
  memory "12 GB"
  //fix input lbael*
  input:
  tuple val(input_label), path(vcf)
   
  output:
  tuple path("BQC19_all.snparray.hwe_filtered.pruned.bed"), path("BQC19_all.snparray.hwe_filtered.pruned.bim"), path("BQC19_all.snparray.hwe_filtered.pruned.fam")
  
  publishDir "${params.result}/bed"

  script:
  """
  ${params.plink2_exec} --vcf ${vcf} --make-bed --out BQC19_all.snparray.hwe_filtered.pruned
  """

}

// step 1 processes each phenotype individually and by type (binary or quantitative) 
process SAIGE_step1 {
  label "SAIGE"
  cache "lenient"
  cpus 3
  scratch '$SLURM_TMPDIR'
  memory "12 GB"

  input:
  tuple val(phen), val(type), path(input_bed), path(input_bim), path(input_fam), path(phen_file)

  output:
  tuple val(phen), val(type), path("${phen}_${type}.rda"), path("${phen}_${type}.varianceRatio.txt")

  publishDir "${params.result}/step1"
  // --skipVarianceRatioEstimation=FALSE \
  script:
  mode = type
  if (mode == 'binary')
  """
      step1_fitNULLGLMM.R \
        --bedFile=${input_bed} \
        --bimFile=${input_bim} \
        --famFile=${input_fam} \
        --phenoFile=${phen_file} \
        --phenoCol=${phen} \
        --sexCol=sex \
        --FemaleCode=2 \
        --MaleCode=1 \
        --covarColList=${params.covar_cols} \
        --qCovarColList=${params.cat_covar} \
        --sampleIDColinphenoFile=BQCID \
        --traitType=${type} \
        --outputPrefix=${phen}_${type} \
        --nThreads=3	\
        --IsOverwriteVarianceRatioFile=TRUE
 
  """
  else
  """
      step1_fitNULLGLMM.R \
        --bedFile=${input_bed} \
        --bimFile=${input_bim} \
        --famFile=${input_fam} \
        --useSparseGRMtoFitNULL=FALSE \
        --phenoFile=${phen_file} \
        --phenoCol=${phen} \
        --sexCol=sex \
        --covarColList=${params.covar_cols} \
        --qCovarColList=${params.cat_covar} \
        --sampleIDColinphenoFile=BQCID \
        --invNormalize=TRUE \
        --traitType=quantitative \
        --outputPrefix=${phen}_${type} \
        --nThreads=3 \
        --IsOverwriteVarianceRatioFile=TRUE
  """

}

//step 2 processes each phen and each chrom separately 
process SAIGE_step2 {
  label "SAIGE"
  cache "lenient"
  cpus 1
  errorStrategy 'ignore'
  memory "10 GB"

  input:
  each chrom
  tuple val(phen), val(type), path(rda), path(varRatio), val(var_type), path(vcf), path(vcf_csi), path(sample_list), path(bqc_males)

  output:
  path("${phen}.${chrom}.${var_type}.SAIGE.gwas.txt")

  publishDir "${params.result}/step2"
  // --IsOutputAFinCaseCtrl=TRUE \
  //--is_output_moreDetails=TRUE
  // --SampleFile=${sample_list} \
  // // --isDropMissingDosages=TRUE
  //         --X_PARregion="10001-2781479,155701383-156030895" \
  //       --is_rewrite_XnonPAR_forMales=TRUE \
  //       --sampleFile_male=${bqc_males} \
  //         --maxMissing=0.01 \

  script:
  mode = type
  if (mode == 'binary')
  """
  step2_SPAtests.R        \
        --vcfFile=${vcf}  \
        --vcfFileIndex=${vcf_csi} \
        --chrom=${chrom} \
        --vcfField=GT \
        --minMAF=0 \
        --minMAC=20 \
        --SAIGEOutputFile=${phen}.${chrom}.${var_type}.SAIGE.gwas.txt \
        --GMMATmodelFile=${rda} \
        --varianceRatioFile=${varRatio} \
        --is_Firth_beta=TRUE    \
        --pCutoffforFirth=0.05 \
        --is_output_moreDetails=TRUE    \
        --sampleFile=${sample_list} \
        --LOCO=TRUE

  """
  else
  """
    step2_SPAtests.R        \
        --vcfFile=${vcf}  \
        --vcfFileIndex=${vcf_csi} \
        --chrom=${chrom} \
        --vcfField=GT \
        --minMAF=0 \
        --minMAC=20 \
        --SAIGEOutputFile=${phen}.${chrom}.${var_type}.SAIGE.gwas.txt \
        --GMMATmodelFile=${rda} \
        --varianceRatioFile=${varRatio} \
        --sampleFile=${sample_list} \
        --LOCO=TRUE
  """
} 


workflow {
// snp array vcf 
input_vcf =  Channel.fromPath(params.inputFiles).map{file -> [file.getBaseName(), file]}

// convert snp array vcf to bed
input_bed_bim_fam = PlinkFormatting(input_vcf)

// read each line of the phen file, fields[0]=phen and fields[1]=type (binary or quant)
phen_type = Channel.from(file(params.inputPhenType).readLines()).map { line -> fields = line.split(); [ fields[0], fields[1] ] }

phen_file = Channel.fromPath(params.inputPhen)

//prepare files for step1 SAIGE
step1_files = input_bed_bim_fam.concat(phen_file).collect()
input_step1 = phen_type.combine(step1_files)

// step1 SAIGE
output_step1 = SAIGE_step1(input_step1)

// for structural variants gives the variant type, vcf, vcf.csi
input_vcf = Channel.fromPath(params.inputvcf).map{vcf_file -> [vcf_file.getName().split('\\.')[0], vcf_file, vcf_file + ".csi" ]}

// list of samples in input_vcf
input_sample_list = Channel.fromPath(params.sampleFile)

// list of male samples 
input_male_sample_list = Channel.fromPath(params.bqc_males)

// prepare files for step2 SAIGE
input_step2 = output_step1.combine(input_vcf)
input_step2_final = input_step2.combine(input_sample_list)
input_step2_final1 = input_step2_final.combine(input_male_sample_list)

//step2 SAIGE
output_step2 = SAIGE_step2(params.chrom, input_step2_final1)
//input_vcf = Channel.fromPath(params.inputvcf).map{vcf_file -> [vcf_file.getName().split('\\.')[0], vcf_file.getName().split('\\.')[3], vcf_file, vcf_file + ".csi" ]}

// add chrom to tuple:
// output_step2 = SAIGE_step2(input_step2_final)

}

