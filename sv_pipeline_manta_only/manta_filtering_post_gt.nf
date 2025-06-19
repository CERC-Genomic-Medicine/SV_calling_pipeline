#!/usr/bin/env nextflow


process MergeGenotyping{
cache "lenient"
cpus 1
memory "50 GB"
//memory "12 GB"
time "20h"
// // maxRetries 2
// // errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
scratch '$SLURM_TMPDIR'
container "${params.smooveContainer}"
containerOptions "-B ${params.referenceDir}:/ref"

input:
tuple val(svtype), path(vcf), path(vcf_csi)

output:
tuple val(svtype), path("${svtype}.smoove.square.vcf.gz")

publishDir "${params.result}/all_batches/merge"
script:

"""
smoove paste --name ${svtype} --outdir . ${vcf} 
"""

}

process DepthFiltering{
cache "lenient"
cpus 1
memory "16 GB"
time "2h"
scratch '$SLURM_TMPDIR'

// bcftools index ${svtype}.duphold.filt.vcf.gz
// +setGT ${vcf} -- -t q -n . -i
input:
tuple val(svtype), path(vcf)

output:
tuple val(svtype), path("${svtype}.duphold.filt.vcf.gz")

publishDir "${params.result}/all_batches/duphold"

script:
mode = svtype
if (mode == 'DEL')
"""
bcftools filter --set-GTs . -e 'GT="alt" & FMT/DHFFC >=0.7' ${vcf} | bcftools view - -Oz -o ${svtype}.duphold.filt.vcf.gz
"""

if (mode == 'DUP')
"""
bcftools filter --set-GTs . -e 'GT="alt" & FMT/DHBFC <= 1.25' ${vcf} | bcftools view - -Oz -o ${svtype}.duphold.filt.vcf.gz
"""

else
"""
bcftools view ${vcf} -Oz -o ${svtype}.duphold.filt.vcf.gz
"""
}

process MissingnessFiltering{
cache "lenient"
cpus 1
memory "16 GB"
time "2h"
scratch '$SLURM_TMPDIR'

input:
tuple val(svtype), path(vcf)

output:
tuple val(svtype), path("${svtype}.missingness.filt.vcf.gz")

publishDir "${params.result}/all_batches/missingness"

script:
mode = svtype
if( mode == 'DEL' || mode == 'DUP' || mode == 'INV') 
"""
bcftools view -i 'F_MISSING < 0.1' ${vcf} | bcftools view - -Oz -o ${svtype}.missingness.filt.vcf.gz
"""
else 
"""
export R_LIBS_USER="/home/oliviagc/projects/rrg-vmooser/oliviagc/R/library"
bcftools +missing2ref ${vcf}| bcftools view - -Oz -o samples_merged_${svtype}.filt2.vcf.gz
bcftools view -i 'F_MISSING < 0.1' samples_merged_${svtype}.filt2.vcf.gz | bcftools view - -Oz -o ${svtype}.missingness.filt.vcf.gz
"""
}



process FindOutliers{
  cache "lenient"
  cpus 1
  memory "80 GB"
  time "20h"
  scratch '$SLURM_TMPDIR'

  input:   
  tuple val(svtype), path(vcf)
  

  output:
  path("outlier_samples.${svtype}*blacklist")
   
  publishDir "${params.result}/all_batches/outliers"
  
  
  script:
  """
  export R_LIBS_USER="/home/oliviagc/projects/rrg-vmooser/oliviagc/R/library"
  # this R script is for one SV type at a time
  Rscript ${params.scripts}/removeSampleOutliers.R ${vcf} samples_merged_${svtype}.filt1.vcf.gz outlier_samples.${svtype}.SM.blacklist

  """
}

process MergeOutliers{
cache "lenient"
cpus 1
memory "32 GB"
time "6h"
scratch '$SLURM_TMPDIR'
input:

path(blacklists)

output:
path("samplesToMerge.list")

publishDir "${params.result}/all_batches/outliers_ids"

script:
"""
cat ${blacklists} | sort | uniq > outlier_samples.blacklist
cat outlier_samples.blacklist | cut -f1 > outlier_ids
bcftools query -l /home/oliviagc/scratch/manta_filtering/results/all_batches/merge/DEL.smoove.square.vcf.gz | grep -v -x -f outlier_ids > samplesToMerge.list
"""
}



process Filter{
cache "lenient"
cpus 1
memory "80 GB"
time "12h"
// container "${params.smooveContainer}"
// containerOptions "-B /home/oliviagc/scratch:/ref"

input:
tuple val(svtype), path(vcf)
path(inclusion_list)

output:
path("samples_merged_${svtype}.Final.vcf.gz")

publishDir "${params.result}/all_batches/final_filtering"

"""
export R_LIBS_USER="/home/oliviagc/projects/rrg-vmooser/oliviagc/R/library"

Rscript ${params.scripts}/renameSVIds.R ${vcf} ${svtype}.filt.vcf.gz
bcftools view --force-samples -S ${inclusion_list} ${svtype}.filt.vcf.gz -Oz > samples_merged_${svtype}.filt2.vcf.gz

# c1 excludes 0/0 genotypes, sort sorts by chromosome number 
bcftools view -c1 samples_merged_${svtype}.filt2.vcf.gz | bcftools sort -Oz -o samples_merged_${svtype}.filt3.vcf.gz 

Rscript ${params.scripts}/renameSVIds.R samples_merged_${svtype}.filt3.vcf.gz samples_merged_${svtype}.filt5.vcf.gz

zcat samples_merged_${svtype}.filt5.vcf.gz | bcftools annotate -x "INFO/SUPP_VEC,INFO/SUPP" > samples_merged_${svtype}.filt6.vcf

java -jar \$EBROOTPICARD/picard.jar FixVcfHeader I=samples_merged_${svtype}.filt6.vcf O=samples_merged_${svtype}.filt7.vcf

bgzip -f samples_merged_${svtype}.filt7.vcf

bcftools +fill-tags -Oz -o samples_merged_${svtype}.Final.vcf.gz samples_merged_${svtype}.filt7.vcf.gz

"""

}


process Annotation {
  cache "lenient"
  cpus 1
  memory "80 GB"
  time "2h"

  input:
  tuple val(sv_type), path(vcf)

  output:
  path("*.Final.clean.annotated.tsv")

  publishDir "${params.result}/all_batches/Annotation"

  // I had to use bed files and cut most fields because AnnotSV memort with TCL (?) is limited 
  script:
  """
  export ANNOTSV="/home/oliviagc/projects/rrg-vmooser/oliviagc/AnnotSV"
  
  gunzip -f ${vcf}
  ${params.survivor} vcftobed ${sv_type}.Final.vcf -99999999 99999999 ${sv_type}.Final.bed
  cat ${sv_type}.Final.bed | cut -f1,2,5,7,11 | awk -F '\t' '{ \$3 = (\$3 == "0" ? \$2+1 : \$3) } 1' OFS='\t'| awk '(\$3 > \$2 )' > ${sv_type}.Final.clean.bed
  /home/oliviagc/projects/rrg-vmooser/oliviagc/AnnotSV/bin/AnnotSV -SVinputFile ${sv_type}.Final.clean.bed -outputDir ./ -SVinputInfo 1 -reciprocal 1 -svtBEDcol 5

  """
}

workflow{
// inputMerge = Channel.fromPath(params.inputFilesMergeGt).map{vcf_file -> [vcf_file.getName().split('\\.')[0],vcf_file.getName().split('\\.')[1], vcf_file,vcf_file + ".csi" ]}
// outputDepthFilter = DepthFiltering(inputDepthFilter).map{vcf_file, csi_file-> [vcf_file.getName().split('\\.')[0], vcf_file, vcf_file + ".csi" ]}.groupTuple(by:0)
// inputMerge = Channel.fromPath(params.inputFilesMergeGt).map{vcf_file -> [vcf_file.getName().split('\\.')[0], vcf_file, vcf_file + ".csi" ]}.groupTuple(by:0)

// //outputMerge = MergeGenotyping(inputMerge) // results in sv_type and vcf.gz 
// outputMerge = Channel.fromPath(params.outputMergeVarRemoved).map{vcf_file -> [vcf_file.getName().split('\\.')[0], vcf_file]}.groupTuple(by:0)
// outputDepthFilter = DepthFiltering(outputMerge)
// inputInsertions = Channel.fromPath(params.inputFilesMantaIns).map{vcf_file -> [vcf_file.getName().split('\\.')[1], vcf_file]}

// outputMissingnessFilter = MissingnessFiltering(outputDepthFilter.concat(inputInsertions))
// outputOutliers = FindOutliers(outputMissingnessFilter).collect()
// inclusionList = MergeOutliers(outputOutliers)
inputMerge = Channel.fromPath(params.outputMissingnessFilter).map{vcf_file -> [vcf_file.getName().split('\\.')[0], vcf_file]}.groupTuple(by:0)
outputFiltering = Filter(inputMerge, params.inclusionList).map{vcf_file -> [vcf_file.getName().split('\\.')[0], vcf_file]}.groupTuple(by:0)
// outputAnnotation = Annotation(outputFiltering)
// MergeAnnotation(outputAnnotation)
}