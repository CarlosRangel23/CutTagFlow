nextflow.enable.dsl = 2

process DIFFBIND {
    tag "DiffBind DBA (${label} - ${histone_mark})"
    label 'process_medium'

    publishDir "${params.outdir}/04_peak_calling/consensus/diffbind/${label}", mode: 'copy'

    input:
    tuple val(meta), path(peaks), path(bams), val(min_overlap)

    output:
    path "*.bed"        , emit: bed
    path "*.RData"      , emit: rdata

    script:
    label        = meta[0]
    histone_mark = meta[1]

    def peaks_arg = peaks.join(',')
    def bams_arg  = bams.join(',')

    """
    diffbind_extraction.R \\
        "${label}" \\
        "${histone_mark}" \\
        "${peaks_arg}" \\
        "${bams_arg}" \\
        ${min_overlap}
    """
}