nextflow.enable.dsl = 2

process BCFTOOLS {
    tag "${histone_mark}"
    label 'process_high'
    publishDir "${params.outdir}/06_variant_calling/${histone_mark}/bcftools", mode: 'copy'

    input:
    tuple val(histone_mark), path(bams), path(bais)
    path fasta
    path fai

    output:
    tuple val(histone_mark),
          path("cuttag_${histone_mark}_bcftools.filtered.vcf.gz"),
          path("cuttag_${histone_mark}_bcftools.filtered.vcf.gz.tbi"), emit: vcf
    tuple val(histone_mark),
          path("cuttag_${histone_mark}_bcftools.stats.txt"), emit: stats
    tuple val(histone_mark),
          path("cuttag_${histone_mark}_bcftools.filtered.bed"),
          path("cuttag_${histone_mark}_bcftools.filtered.bim"),
          path("cuttag_${histone_mark}_bcftools.filtered.fam"), emit: plink

    script:
    def prefix = "cuttag_${histone_mark}_bcftools"
    """
    set -o pipefail

    bcftools mpileup \
        -f ${fasta} \
        -q 20 \
        -Q 20 \
        -d 250 \
        --no-BAQ \
        ${bams} \
        -Ou | \
    bcftools call \
        --threads ${task.cpus} \
        -mv \
        -Oz \
        -o ${prefix}.raw.vcf.gz

    bcftools filter \
        -e 'QUAL<20 || DP<10 || MQ<30' \
        ${prefix}.raw.vcf.gz \
        -Oz \
        -o ${prefix}.filtered.vcf.gz

    tabix -p vcf ${prefix}.filtered.vcf.gz

    bcftools stats ${prefix}.filtered.vcf.gz > ${prefix}.stats.txt

    plink2 \
        --vcf ${prefix}.filtered.vcf.gz \
        --make-bed \
        --out ${prefix}.filtered \
        --autosome \
        --snps-only \
        --max-alleles 2
    """
}
