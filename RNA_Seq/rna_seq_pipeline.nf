#!/usr/bin/env nextflow

process Fastp {
  cache "lenient"
  cpus 3
  scratch '$SLURM_TMPDIR'
  maxRetries 2
  errorStrategy "ignore"
  memory "12 GB"
  time "6h"
  
  input:
  tuple val(var_id), val(bqc_id), path(R1), path(R2)

  output:
  tuple path("*R1.fq.gz"), path("*R2.fq.gz"), path("fastp.html"), val(bqc_id), val(var_id)
  
  publishDir "${params.result}/rna_seq_fastp/${bqc_id}/${var_id}", pattern: "*.R*.fq.gz", mode: 'copy'

  """
  module load fastp
  fastp --qualified_quality_phred 20 --unqualified_percent_limit 20 --correction -i ${R1} -I ${R2} -o out_fastp_${bqc_id}.${var_id}.R1.fq.gz -O out_fastp_${bqc_id}.${var_id}.R2.fq.gz
  """


}

process FastQC {
  cache "lenient"
  cpus 5
  scratch '$SLURM_TMPDIR'
  memory "40 GB"
  time "6h"
  
  input:
  tuple path(R1), path(R2), path(html), val(bqc_id), val(var_id)


  output:
  path("*fastqc*")
  
  publishDir "${params.result}/rna_seq_fastqc_step2/${bqc_id}/${var_id}", pattern:"*fastqc*", mode: "copy"

  """
  module load fastqc
  
  fastqc ${R1} ${R2}
  
  """


}



process Mapping_hg38 {
  cache "lenient"
  cpus 5
  scratch '$SLURM_TMPDIR'
  errorStrategy "ignore"
  memory "35 GB"
  time "6h"
  
  input:
  tuple val(var_id), val(bqc_id), path(R1), path(R2)

  output:
  tuple path(R1), path(R2), path("*.bam"), path("*.final.out")
  
  publishDir "${params.result}/rna_seq_star_step3/${bqc_id}/${var_id}", pattern: "*.bam", mode: "copy"
  publishDir "${params.result}/rna_seq_star_step3/${bqc_id}/${var_id}/summary", pattern: "*.final.out", mode: "copy"

  script:
  """
  module load star
  STAR --runThreadN 5 \
  --outFilterMatchNmin 90 \
  --outFilterMatchNminOverLread 0.95 \
  --genomeDir ~/scratch/STAR_genomic_indices_gencode \
  --readFilesIn ${R1}, ${R2} \
  --outFileNamePrefix ${bqc_id}_${var_id} \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM SortedByCoordinate
  
  """


}

process Annotation_hg38 {
  cache "lenient"
  cpus 3
  // maxRetries 2
  errorStrategy "ignore"
  memory "16 GB"
  time "12h"
  
  input:
  tuple val(bqc), val(var), path(bam_file)

  output:
  path("${bqc}_${var}.txt*")
  
  publishDir "${params.result}/annotation_step4_rerun/", pattern:"${bqc}_${var}.txt*", mode:"copy"

  """
  module load subread 
  featureCounts -T 3 -p --countReadPairs \
  -t exon \
  -g gene_id \
  -a ~/scratch/gencode.v47.annotation.gtf \
  -o ${bqc}_${var}.txt ${bam_file}
  
  """


}

process DESeq2 {
  cache "lenient"
  cpus 5
  memory "40 GB"
  time "6h"
  
  input:
  path(txt_file_all_samples)

  output:
  path("vst_mat.csv")

  publishDir "${params.result}/normalization_step5", mode: "copy"
  script:
  """
  module load r 
  export R_LIBS_USER="/home/oliviagc/projects/rrg-vmooser/oliviagc/R/library"
  Rscript /home/oliviagc/scratch/RNA_Seq/deseq.R ${txt_file_all_samples}

  """
}


workflow {

input_rna_seq_fastp_r1 =  Channel.fromPath(params.input_rna_seq_fastp_r1)
                                  .map{file -> [file.getName().split('\\.')[1].split('\\-')[0],file.getName().split('\\.')[0], file]}
input_rna_seq_fastp_r2 =  Channel.fromPath(params.input_rna_seq_fastp_r2)
                                  .map{file -> [file.getName().split('\\.')[1].split('\\-')[0], file.getName().split('\\.')[0], file]}

fastp_output=  input_rna_seq_fastp_r1.join(input_rna_seq_fastp_r2, by: [ 0, 1])

// star_mapping_output_hg38 = Mapping_hg38(fastp_output)
star_mapping_output_hg38 = Channel.fromPath(params.output)
                                  .map{file -> [file.getName().split('\\_')[2],file.getName().split('\\_')[3], file]}

annotation_output_hg38 = Annotation_hg38(star_mapping_output_hg38)

// annotation_output_hg38 = Channel.fromPath(params.merged_output)
// normalized = DESeq2(annotation_output_hg38)



}