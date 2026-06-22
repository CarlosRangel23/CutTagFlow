nextflow.enable.dsl = 2

process FILTER_BAM {

    tag "${histone_mark}/${sample} (${label})"
    label 'process_medium'

    input:
    tuple val(sample), path(bam), path(bai), val(histone_mark)
    val(label)
    path blacklist_bed

    output:
    tuple val(sample), path("${sample}.${label}.filtered.bam"), val(label), val(histone_mark), emit: bam
    path "${sample}.${label}.filtered.bam.bai",                                                emit: bai
    path "${sample}.${label}.idxstats.txt",                                                    emit: idxstats
    path "${sample}.${label}.fragmentLen.txt",                                                 emit: frag_len

    publishDir "${params.outdir}/03_filtered/${histone_mark}/${label}/${sample}", mode: 'copy'

    script:
    def samtools_flags = (label == 'noDups') ? "-F 1804" : "-F 780"
    def chrs = "chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

    """
    samtools view -h -q 2 ${samtools_flags} ${bam} $chrs | \\
    samtools sort -O bam - | \\
    bedtools intersect -v -abam stdin -b ${blacklist_bed} > ${sample}.${label}.filtered.bam

    samtools index ${sample}.${label}.filtered.bam
    samtools idxstats ${sample}.${label}.filtered.bam > ${sample}.${label}.idxstats.txt

    samtools view -F 0x04 ${sample}.${label}.filtered.bam | \\
      awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs(\$9)}' | \\
      sort | uniq -c | \\
      awk -v OFS="\\\\t" '{print \$2, \$1/2}' > ${sample}.${label}.fragmentLen.txt
    """
}