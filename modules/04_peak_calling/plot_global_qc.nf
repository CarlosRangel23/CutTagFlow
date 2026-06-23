nextflow.enable.dsl = 2

process PLOT_GLOBAL_QC {
    tag "Global QC Plotting"
    label 'process_low'
    publishDir "${params.outdir}/05_global_qc", mode: 'copy'

    input:
    path qc_files       // Receives ALL .frip.csv and .tsse.csv files
    path cutoff_files   // Receives ALL _cutoff_analysis.txt files from MACS3

    output:
    path "FRiP_peaks_comparison.png"     , emit: frip_peaks_plot
    path "FRiP_TSS_comparison.png"       , emit: frip_tss_plot
    path "TSSE_scores_comparison.png"    , emit: tsse_plot
    path "MACS3_cutoff_peaks_impact.png" , emit: cutoff_plot
    path "Experiment_QC_Summary.csv"     , emit: final_table

    publishDir "${params.outdir}/04_peak_calling/Advanced_QC", mode: 'copy'

    script:
    """
    plot_advanced.qc.R
    """
}