nextflow.enable.dsl = 2

process MARK_DUPLICATES {

    tag "$sample"
    label 'process_medium'

    input:
    tuple val(sample), path(bam), val(histone_mark)

    output:
    tuple val(sample), path("${sample}.sorted.dupMarked.bam"), val(histone_mark), emit: bam
    path "${sample}.sorted.dupMarked.bam.bai"                                   , emit: bai
    path "${sample}_picard.rmDup.txt"                                           , emit: picard_metrics
    path "${sample}_fragmentLen.txt"                                            , emit: frag_len
    path "${sample}_idxstats.txt"                                               , emit: idxstats

    publishDir "${params.outdir}/02_alignment/temp_picard_idxstats/${sample}", mode: 'copy'

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