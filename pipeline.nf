#!/usr/bin/env nextflow 

process DellySV {
// call Delly, compress and index resulting file 
   cache "lenient"
   cpus 6
   memory "20 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

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
   scratch '$SLURM_TMPDIR'

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.manta.vcf.gz"), path("${input_label}.manta.vcf.gz.tbi")

   publishDir "${params.result}/manta", pattern: '*.vcf.gz*', mode: 'copy'

   """
   # call manta, restricted to ref regions included in bed file 
   ${params.manta}/bin/configManta.py --bam ${cram} --referenceFasta ${params.referenceGenome} --callRegions ${params.referenceDir}/GRCh38.include.bed.gz
   # execute script with specified job count equal to 4 
   ./MantaWorkflow/runWorkflow.py -j 4
   
   # analysis of BNDs to determine INVs, output = vcf
   ${params.manta}/libexec/convertInversion.py ${params.manta}/libexec/samtools ${params.referenceGenome} ./MantaWorkflow/results/variants/diploidSV.vcf.gz > ${input_label}.manta.vcf
   
   #compression to vcf.gz 
   bcftools view ${input_label}.manta.vcf -Oz -o ${input_label}.manta.vcf.gz
   
   #tabix index
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
   samtools view -h -b -F 1294 ${cram} chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y} > ${input_label}.discordant.unsorted.bam
   
   samtools view -h ${cram} chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y} | ${params.lumpy}/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${input_label}.splitters.unsorted.bam
   
   samtools sort -o ${input_label}.discordant.bam ${input_label}.discordant.unsorted.bam
   samtools sort -o ${input_label}.splitters.bam ${input_label}.splitters.unsorted.bam 

   sed "s|SAMTOOLS=|SAMTOOLS=\$(which samtools)|" ${params.lumpy}/bin/lumpyexpress.config > lumpyexpress.config
   sed -i "s|PYTHON=.*|PYTHON=\$(which python)|" lumpyexpress.config
   
   ${params.lumpy}/bin/lumpyexpress -B ${cram} -S ${input_label}.discordant.bam -D ${input_label}.splitters.bam -o ${input_label}.vcf -K lumpyexpress.config 
   cat ${input_label}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > ${input_label}.lumpy.vcf
   bgzip ${input_label}.lumpy.vcf
   tabix -p vcf ${input_label}.lumpy.vcf.gz
   """
 }

 
 process CNVnator{

   cache "lenient"
   cpus 4
   memory "16 GB"
   time "12h"
   scratch '$SLURM_TMPDIR'

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.cnvnator.vcf.gz"), path("${input_label}.cnvnator.vcf.gz.tbi")

   publishDir "${params.result}/cnvnator", pattern: '*.vcf.gz*', mode: 'copy'
   // module load bcftools, samtools, gcc, root 
   script:
   """
   #create ROOTSYS variable 
   export ROOTSYS="./root"
   export PATH="\$PATH:\$ROOTSYS/bin"
   export LD_LIBRARY_PATH="\${LD_LIBRARY_PATH:-}:\$ROOTSYS/lib"

   ${params.cnvnator}/cnvnator -root ${input_label}.root -tree ${cram} -chrom chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

   ${params.cnvnator}/cnvnator -root ${input_label}.root -his 200 -fasta ${params.referenceGenome} 
   
   ${params.cnvnator}/cnvnator -root ${input_label}.root -stat 200

   ${params.cnvnator}/cnvnator -root ${input_label}.root -partition 200

   ${params.cnvnator}/cnvnator -root ${input_label}.root -call 200 > ${input_label}.cnvnator


   perl ${params.cnvnator}/cnvnator2VCF.pl -reference GRCh37 ${input_label}.cnvnator ${params.referenceGenome} > ${input_label}.cnvnator.vcf


   # compress
   bcftools view ${input_label}.cnvnator.vcf -Oz -o ${input_label}.cnvnator.vcf.gz
   
   # indexation
   tabix -p vcf ${input_label}.cnvnator.vcf.gz
   """

 }

process BreakDancer{
   cache "lenient"
   cpus 4
   memory "16 GB"
   time "12h"
   scratch '$SLURM_TMPDIR'
   
   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.ctx"), path("${input_label}.breakdancer.vcf.gz"), path("${input_label}.breakdancer.vcf.gz.tbi")

   publishDir "${params.result}/breakdancer", pattern: '*.vcf.gz*', mode: 'copy'
   publishDir "${params.result}/breakdancer/ctx", pattern: '*.ctx', mode: 'copy'

   script:
   // module load bcftools, samtools, perl, require GD::Graph::histogram, Statistics::Descriptive, Math::CDF perl packages 
   // remove the variants where POS > END from vcf file 
   """
   export PATH="${params.perl}/bin:\$PATH"
   export PERL5LIB="${params.perl}/lib/perl5:\$PERL5LIB"
   export PERL_LOCAL_LIB_ROOT="${params.perl}"
   export PERL_MB_OPT="--install_base "${params.perl}""
   export PERL_MM_OPT="INSTALL_BASE=${params.perl}"
   
   samtools view -h -b -o ${input_label}.bam ${cram} chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}
   ${params.breakdancer}/perl/bam2cfg.pl ${input_label}.bam > ${input_label}.cfg
   ${params.breakdancer}/install/bin/breakdancer-max -q 30 -n 10000 ${input_label}.cfg > ${input_label}.ctx
   python ${params.scripts}/breakdancer2vcf.py ${params.referenceGenome} ${input_label}.ctx ${input_label}.vcf

   # compress
   bcftools filter -e 'POS > INFO/END' ${input_label}.vcf -Oz -o ${input_label}.breakdancer.vcf.gz
   
   # indexation
   bcftools index -t ${input_label}.breakdancer.vcf.gz
   """

}

process MeltALU {
   cache "lenient"
   cpus 3
   memory "12 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.ALU.tar.gz") 

   publishDir "${params.result}/melt/ALU", pattern: '*.tar.gz', mode: 'copy'
   
   script:
   """
   mkdir -p ./ALU
   java -Xmx8G -jar ${params.melt}/MELT.jar Preprocess -bamfile ${cram} -h ${params.referenceGenome}
   java -Xmx12G -jar ${params.melt}/MELT.jar IndivAnalysis -bamfile ${cram} -w ./ALU/ -t ${params.mei}/ALU_MELT.zip -h ${params.referenceGenome}
   tar -zcvf ${input_label}.ALU.tar.gz ALU 
   """
}

process MeltHERVK {
   cache "lenient"
   cpus 3
   memory "12 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.HERVK.tar.gz") 

   publishDir "${params.result}/melt/HERVK", pattern: '*.tar.gz', mode: 'copy'
   
   script:
   """
   mkdir -p ./HERVK
   java -Xmx8G -jar ${params.melt}/MELT.jar Preprocess -bamfile ${cram} -h ${params.referenceGenome}
   java -Xmx12G -jar ${params.melt}/MELT.jar IndivAnalysis -bamfile ${cram} -w ./HERVK/ -t ${params.mei}/HERVK_MELT.zip -h ${params.referenceGenome}
   tar -zcvf ${input_label}.HERVK.tar.gz HERVK 
   """
}

process MeltLINE1 {
   cache "lenient"
   cpus 3
   memory "12 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.LINE1.tar.gz") 

   publishDir "${params.result}/melt/LINE1", pattern: '*.tar.gz', mode: 'copy'
   
   script:
   """
   mkdir -p ./LINE1
   java -Xmx8G -jar ${params.melt}/MELT.jar Preprocess -bamfile ${cram} -h ${params.referenceGenome}
   java -Xmx12G -jar ${params.melt}/MELT.jar IndivAnalysis -bamfile ${cram} -w ./LINE1/ -t ${params.mei}/LINE1_MELT.zip -h ${params.referenceGenome}
   tar -zcvf ${input_label}.LINE1.tar.gz LINE1 
   """
}

process MeltSVA {
   cache "lenient"
   cpus 3
   memory "12 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:
   tuple val(input_label), file(cram), file(crai)

   output:
   tuple val(input_label), path("${input_label}.SVA.tar.gz") 

   publishDir "${params.result}/melt/SVA", pattern: '*.tar.gz', mode: 'copy'
   
   script:
   """
   mkdir -p ./SVA
   java -Xmx8G -jar ${params.melt}/MELT.jar Preprocess -bamfile ${cram} -h ${params.referenceGenome}
   java -Xmx12G -jar ${params.melt}/MELT.jar IndivAnalysis -bamfile ${cram} -w ./SVA/ -t ${params.mei}/SVA_MELT.zip -h ${params.referenceGenome}
   tar -zcvf ${input_label}.SVA.tar.gz SVA 
   """
}

workflow {
   
   inputFiles = Channel.fromPath(params.inputFilesSubset).map{file -> [file.getSimpleName(), file, file + ".crai"]}
   // inputFiles = Channel.fromPath(params.inputFilesBam).map{file -> [file.getSimpleName(), file, file + ".bai"]}
                 
   inputFiles | MantaSV 
   inputFiles | DellySV 
   inputFiles | LumpySV 
   inputFiles | CNVnator
   inputFiles | BreakDancer
   inputFiles | MeltALU
   inputFiles | MeltHERVK
   inputFiles | MeltLINE1
   inputFiles | MeltSVA
}