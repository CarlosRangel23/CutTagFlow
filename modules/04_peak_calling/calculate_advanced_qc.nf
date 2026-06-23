nextflow.enable.dsl = 2

process ADVANCED_QC {
    tag "${sample} - ${label}"
    label 'process_medium'

    input:
    tuple val(sample), val(histone_mark), path(bam), path(bai), path(peak)
    val label

    output:
    tuple val(sample), val(histone_mark), val(label), path("${sample}.${label}.qc_metrics.csv"), emit: metrics

    script:
    """
    advanced.qc.R 

    """
}