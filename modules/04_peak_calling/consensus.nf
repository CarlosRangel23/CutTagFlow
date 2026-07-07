nextflow.enable.dsl = 2

process CONSENSUS {
    tag "Consensus Decay (${label} - ${histone_mark})"
    label 'process_medium'

    publishDir "${params.outdir}/04_peak_calling/consensus/png/${label}", mode: 'copy'

    input:
    tuple val(meta), path(peaks), path(bams) // Captures: [ [label, mark], [peaks...], [bams...] ]

    output:
    path "*.png", emit: plots

    script:
    // Extract values cleanly within the script context to satisfy tag/publishDir directives
    label        = meta[0]
    histone_mark = meta[1]

    def peaks_arg = peaks.join(',')
    def bams_arg  = bams.join(',')

    """
    consensus_decay.R \\
        "${label}" \\
        "${histone_mark}" \\
        "${peaks_arg}" \\
        "${bams_arg}"
    """
}