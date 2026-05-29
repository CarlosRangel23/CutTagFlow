nextflow.enable.dsl = 2

process PLOT_QC_METRICS {
    tag "ggplot_qc_summary"
    label 'process_low'

    publishDir "${params.outdir}/02_alignment/alignment_summary", mode: 'copy'

    input:
    path picard_files   
    path frag_files     
    path idxstats_files 

    output:
    path "alignment_dup_qc_summary.txt"                 , emit: report
    path "dup_plots.pdf"                                , emit: dup_pdf
    path "mito_sex_plots.pdf"                           , emit: mito_sex_pdf
    path "fragmentlengthplots.pdf"                      , emit: frag_pdf
    path "fragment_length_highlight_all_samples.pdf"    , emit: frag_highlighted_pdf

    script:
    """
    plot_qc_metrics.R
    """
}