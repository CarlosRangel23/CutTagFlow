nextflow.enable.dsl = 2

process GATK_JOINT_GENOTYPING {
    tag "${histone_mark} samples"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/${histone_mark}", mode: 'copy'

    input:
    tuple val(histone_mark), path(gvcfs), path(tbis)

    path fasta
    path fai
    path dict

    output:
    tuple val(histone_mark), path("cuttag_${histone_mark}_joint_genotyped.vcf.gz")    , emit: vcf
    tuple val(histone_mark), path("cuttag_${histone_mark}_joint_genotyped.vcf.gz.tbi"), emit: vcf_index

    script:
    """
    for gvcf in ${gvcfs}; do
        echo "-V \$gvcf" >> gvcfs_flags.list
    done

    gatk --java-options "-Xmx16g" CombineGVCFs \
        -R ${fasta} \
        \$(cat gvcfs_flags.list) \
        -O "${histone_mark}_combined.g.vcf.gz"

    gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R ${fasta} \
        -V "${histone_mark}_combined.g.vcf.gz" \
        -O "cuttag_${histone_mark}_joint_genotyped.vcf.gz"
    """
}