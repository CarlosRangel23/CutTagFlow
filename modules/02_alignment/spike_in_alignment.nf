nextflow.enable.dsl = 2

process ALIGNMENT {
  
    tag "$sample"
    label 'process_high' 

    input:
    tuple val(sample), path(fastq1), path(fastq2)
    path index_dir

    output:
    tuple val(sample), path("${sample}.sorted.bam"), emit: bam
    path "${sample}.sorted.bam.bai"         , emit: bai
    path "${sample}.bowtie2.txt"            , emit: summary

    publishDir "${params.outdir}/02_alignment/${sample}", mode: 'copy'
    publishDir "${params.outdir}/02_alignment/${sample}", mode: 'copy', 

    script:
    def cores = task.cpus
    def ref = "${index_dir}/genome"

    """
    bowtie2 \
      --end-to-end \
      --very-sensitive \
      --no-mixed \
      --no-discordant \
      --phred33 \
      --dovetail \
      -I 10 \
      -X 700 \
      -p ${cores} \
      -x ${ref} \
      -1 ${fastq1} \
      -2 ${fastq2} \
    | samtools sort -@ ${cores} -o ${sample}.sorted.bam -

    samtools index ${sample}.sorted.bam

    cp .command.log ${sample}.bowtie2.txt || true
    
    echo "[ALIGNMENT] Aligned bam generated for ${sample}"
    """
}