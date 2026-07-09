nextflow.enable.dsl = 2

process GATK_JOINT_GENOTYPING {
    tag "${histone_mark} samples"
    label 'process_high'

    publishDir "${params.outdir}/06_variant_calling/cohort_vcf/${histone_mark}", mode: 'copy'

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
    # 1. Generamos la lista de flags -V para los archivos de esta marca
    for gvcf in ${gvcfs}; do
        echo "-V \$gvcf" >> gvcfs_flags.list
    done

    # 2. Combinamos los GVCFs específicos de esta marca de histona
    gatk --java-options "-Xmx16g" CombineGVCFs \
        -R ${fasta} \
        \$(cat gvcfs_flags.list) \
        -O "${histone_mark}_combined.g.vcf.gz"

    # 3. Genotipado conjunto final de la marca
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
        -R ${fasta} \
        -V "${histone_mark}_combined.g.vcf.gz" \
        -O "cuttag_${histone_mark}_joint_genotyped.vcf.gz"
    """
}