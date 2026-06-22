nextflow.enable.dsl = 2

process PLOT_QC_FILTERED {

    label 'process_low'

    input:
    path frag_files       
    path idxstats_files   
    val label

    output:
    path "mito_sex_plots.pdf"                        , emit: mito_sex_pdf
    path "fragmentlengthplots.pdf"                   , emit: frag_len_pdf
    path "fragment_length_highlight_all_samples.pdf" , emit: highlight_pdf

    publishDir "${params.outdir}/03_filtered/global_plots/${label}", mode: 'copy'

    script:
    """
    plot_filtered_qc_metrics.R --label ${label} \
    --frag ${frag_files.join(',')} \
    --idx ${idxstats_files.join(',')}
    """
}