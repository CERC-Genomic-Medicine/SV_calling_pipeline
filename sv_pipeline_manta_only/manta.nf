#!/usr/bin/env nextflow
// Rscripts use vcfR library, which needs r/4.4.0 (StdEnv/2023)
// smoove needs (*** must visit the correct order of loading):
//      - python packages/tools:
//          -pysam (0.15.0 or later), numpy (1.10.2), scipy (0.17.0)
//          -svtyper, which needs scipy and py2.7.x 
//          -all must be compatible versions and work with py2.7.18 
//      - bcftools, samtools, gcc 
//      - lumpy, lumpy filter (add to $PATH)
//      - gsort (import using go)


process FilterSV {
  // for each SVtool, filter the vcf file to remove any SVs smaller than 50 bp and any SVs without 'PASS'
  // separate the vcf file into different files by SVTYPE
  cache "lenient"
  cpus 1
  memory "8 GB"
  time "4h"
  maxRetries 2
  errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
  scratch '$SLURM_TMPDIR'
  
  input:
  tuple val(input_label), val(sv_tool), file(vcf), file(vcf_tbi)
   

  output:
   // For each sv_tool, output one vcf for each sv_type 
  path("${input_label}.*.*.vcf")

  publishDir "${params.result}/filter"

  script:
  """
  export R_LIBS_USER="/home/oliviagc/projects/rrg-vmooser/oliviagc/R/library"
  variants=("BND" "DEL" "DUP" "INS" "INV")

  bcftools view ${vcf} -r ${params.chrom} > ${input_label}.manta01.vcf      
  ${params.survivor} filter ${input_label}.manta01.vcf NA 50 -1 0 -1 ${input_label}.manta02.vcf
  bcftools view -f PASS  ${input_label}.manta02.vcf -Ov -o ${input_label}.manta03.vcf
  Rscript ${params.scripts}/renameSVIds.R ${input_label}.manta03.vcf ${input_label}.manta.ALL.vcf.gz
  gunzip ${input_label}.manta.ALL.vcf.gz

  for val in \${variants[@]}; 
  do
    bcftools view -i 'INFO/SVTYPE="'\${val}'"' ${input_label}.manta.ALL.vcf > ${input_label}.manta.\${val}.vcf
  done

  sed 's/SVTYPE=BND/SVTYPE=TRA/' ${input_label}.manta.ALL.vcf | bcftools view -i 'INFO/SVTYPE="TRA"' > ${input_label}.manta.BNDtoTRA.vcf
  """

}


process MergeSamples{
cache "lenient"
cpus 5
memory "32 GB"
time "4h"
maxRetries 1
errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
scratch '$SLURM_TMPDIR'

input:
tuple val(svtype), path(vcf)

output:
tuple val(svtype), path("samples_merged.${svtype}.fixed.vcf"), path("samples_merged.${svtype}.sites.vcf")
publishDir "${params.result}/merge"

script:

"""
export R_LIBS_USER="/home/oliviagc/projects/rrg-vmooser/oliviagc/R/library"
ls ${vcf} > vcfCallFiles_${svtype}.list
${params.survivor} merge vcfCallFiles_${svtype}.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} callsMerged_${svtype}.vcf
bgzip callsMerged_${svtype}.vcf
bcftools filter -e 'POS > INFO/END' callsMerged_${svtype}.vcf.gz -Oz -o callsMerged_${svtype}_filt.vcf.gz

java -Xmx31g -jar \$EBROOTPICARD/picard.jar SortVcf I=callsMerged_${svtype}_filt.vcf.gz O=callsMerged_${svtype}_sorted.vcf.gz MAX_RECORDS_IN_RAM=10000

# clean up FORMAT section for genotyping
bcftools annotate -x FORMAT/PSV,FORMAT/LN,FORMAT/DR,FORMAT/ST,FORMAT/QV,FORMAT/TY,FORMAT/ID,FORMAT/RAL,FORMAT/AAL,FORMAT/CO callsMerged_${svtype}_sorted.vcf.gz > callsMerged_${svtype}_sorted_clean.vcf
# add rename SV_IDs ?

${params.survivor} filter callsMerged_${svtype}_sorted_clean.vcf NA 50 50000000 0 -1 samples_merged.${svtype}.filt.vcf
Rscript ${params.scripts}/fixCIPOS2.R samples_merged.${svtype}.filt.vcf samples_merged.${svtype}.fixed.vcf.gz
gunzip samples_merged.${svtype}.fixed.vcf.gz
grep "^##" samples_merged.${svtype}.fixed.vcf > header.${svtype}.txt
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> header.${svtype}.txt
bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\n" samples_merged.${svtype}.fixed.vcf  >> header.${svtype}.txt
bgzip header.${svtype}.txt
tabix -p vcf header.${svtype}.txt.gz
bcftools view header.${svtype}.txt.gz -r ${params.chrom}  > samples_merged.${svtype}.sites.vcf


"""

}

process Genotyping{
cache "lenient"
cpus 1
memory "8 GB"
time "12h"
maxRetries 2
errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }
scratch '$SLURM_TMPDIR'

container "${params.smooveContainer}"
containerOptions "-B ${params.referenceDir}:/ref"

input:
tuple val(input_label), path(cram), path(crai), val(svtype), path(vcf_sites)

output:
tuple val(svtype), path("*smoove.genotyped.vcf*")

publishDir "${params.result}/genotype"
script:

"""
smoove genotype -p 1 --name ${svtype}.${input_label} --outdir . --fasta /ref/Homo_sapiens.GRCh38.fa --duphold --vcf ${vcf_sites} ${cram}
"""

}

process MergeGenotyping{
cache "lenient"
cpus 1
memory "4 GB"
time "12h"
maxRetries 2
errorStrategy { task.attempt <= maxRetries ? 'retry' : 'ignore' }

container "${params.smooveContainer}"
containerOptions "-B ${params.referenceDir}:/ref"

input:
tuple val(svtype), path(vcf), path(vcf_csi)

output:
tuple val(svtype), path("*")

publishDir "${params.result}/mergegenotype"
script:

"""
smoove paste --name ${svtype} --outdir . ${vcf}
"""

}

workflow {

Channel.of(["DEL"],["DUP"],["INV"]) \
 | set { to_keep }
inputFilesFilter =  Channel.fromPath(params.inputFilesFilter).map{file -> [file.getName().split('\\.')[0], file.getName().split('\\.')[1], file, file + ".tbi"]}

//inputMerge = FilterSV(inputFilesFilter).flatMap().map{vcf_file -> [vcf_file.getName().split('\\.')[2], vcf_file]}.groupTuple(by:0)

//inputMerge = Channel.fromPath(params.inputMergeFiles).map{vcf_file -> [vcf_file.getName().split('\\.')[2], vcf_file]}.groupTuple(by:0)
//outputMerge = MergeSamples(inputMerge).join(to_keep)

outputMerge = Channel.fromPath(params.outputMergeFiles).map{vcf_file -> [vcf_file.getName().split('\\.')[1], vcf_file]}.groupTuple(by:0).join(to_keep)
inputGenotype = Channel.fromPath(params.inputCram).map{file -> [file.getSimpleName(), file, file + ".crai"]}
outputGenotype = Genotyping(inputGenotype.combine(outputMerge))


// inputMergeGenotype = Channel.fromPath("/home/oliviagc/scratch/mantaAnalysis/results/Batch8Test/genotype/*.*-smoove.genotyped.vcf.gz").map{vcf_file -> [vcf_file.getName().split('\\.')[0], vcf_file,vcf_file + ".csi" ]}.groupTuple(by:0).view()
// outputGenotypeMerge = MergeGenotyping(inputMergeGenotype)
}

