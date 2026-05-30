process REMOVE_DUPLICATES {

    tag "$sample"
    label 'process_medium'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}.dedup.bam"), emit: bam
    path "${sample}.dedup.bam.bai", emit: bai
    path "${sample}_dedup_idxstats.txt", emit: idxstats
    path "${sample}_dedup_fragmentLen.txt", emit: frag_len

    publishDir "${params.outdir}/02_alignment/dedup/${sample}", mode: 'copy'

    script:
    """
    samtools view -b -F 1024 $bam > ${sample}.dedup.bam
    samtools index ${sample}.dedup.bam

    samtools idxstats ${sample}.dedup.bam > ${sample}_dedup_idxstats.txt

    samtools view -F 0x04 ${sample}.dedup.bam | \
      awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs(\$9)}' | \
      sort | uniq -c | \
      awk -v OFS="\\t" '{print \$2, \$1/2}' > ${sample}_dedup_fragmentLen.txt
    """
}