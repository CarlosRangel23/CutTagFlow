nextflow.enable.dsl = 2

process PLOT_ALIGNMENT {
    tag "ggplot_alignment"
    label 'process_low' 

    input:
    path(summary) 

    output:
    path "alignment_efficiency.pdf", emit: pdf

    publishDir "${params.outdir}/02_alignment/plot", mode: 'copy'

    script:
    """
    plotting_alignment.R
    """
}