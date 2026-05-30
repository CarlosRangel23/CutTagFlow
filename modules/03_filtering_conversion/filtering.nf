process FILTER_BAM {

    tag "$sample ($label)"

    input:
    tuple val(sample), path(bam), val(label)

    output:
    tuple val(sample), path("${sample}.${label}.filtered.bam"), val(label), emit: bam
    path "${sample}.${label}.idxstats.txt", emit: idxstats
    path "${sample}.${label}.fragmentLen.txt", emit: frag_len

    publishDir "${params.outdir}/02_alignment/filtered/${label}/${sample}", mode: 'copy'

    script:
    """
    samtools view -b -q 30 -f 2 -F 1804 $bam > ${sample}.${label}.filtered.bam
    samtools index ${sample}.${label}.filtered.bam

    samtools idxstats ${sample}.${label}.filtered.bam > ${sample}.${label}.idxstats.txt

    samtools view -F 0x04 ${sample}.${label}.filtered.bam | \
      awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs(\$9)}' | \
      sort | uniq -c | \
      awk -v OFS="\\t" '{print \$2, \$1/2}' > ${sample}.${label}.fragmentLen.txt
    """
}