nextflow.enable.dsl = 2

process PLOT_ALIGNMENT {

    tag "ggplot_alignment"
    label 'process_low' 
    publishDir "${params.outdir}/02_alignment/alignment_summary", mode: 'copy'

    input:
    path(summary) 

    output:
    path "alignment_summary_plots.pdf", emit: summary_pdf
    path "alignment_summary_report.txt", emit: report


    script:
    """
    plotting_alignment.R ${summary}
    """
}