nextflow.enable.dsl = 2

process MARK_DUPLICATES {

    tag "$sample"
    label 'process_medium'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.sorted.dupMarked.bam"), emit: bam
    path "${sample}.sorted.dupMarked.bam.bai"                , emit: bai
    path "${sample}_picard.rmDup.txt"                        , emit: picard_metrics
    path "${sample}_fragmentLen.txt"                         , emit: frag_len
    path "${sample}_idxstats.txt"                            , emit: idxstats

    publishDir "${params.outdir}/03_filtered/summary_post_alignment/${sample}", mode: 'copy'

    script:
    """
    picard MarkDuplicates \
      I=${bam} \
      O=${sample}.sorted.dupMarked.bam \
      REMOVE_DUPLICATES=false \
      METRICS_FILE=${sample}_picard.rmDup.txt

    samtools index ${sample}.sorted.dupMarked.bam
    samtools idxstats ${sample}.sorted.dupMarked.bam > ${sample}_idxstats.txt
    
    samtools view -F 0x04 ${sample}.sorted.dupMarked.bam | \
      awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs(\$9)}' | \
      sort | uniq -c | \
      awk -v OFS="\t" '{print \$2, \$1/2}' > ${sample}_fragmentLen.txt    
    """
}