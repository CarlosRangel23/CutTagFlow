nextflow.enable.dsl = 2

process GATK_HAPLOTYPE_CALLER {
    tag "${sample}"
    label 'process_medium'

    publishDir "${params.outdir}/06_variant_calling/gvcf", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai), val(histone_mark)

    path fasta
    path fai
    path dict

    output:
    tuple val(histone_mark), path("${sample}.g.vcf.gz"), path("${sample}.g.vcf.gz.tbi"), emit: gvcf

    script:    
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -O "${sample}.g.vcf.gz" \
        -ERC GVCF \
        --sample-name "${sample}" \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation
    """
}