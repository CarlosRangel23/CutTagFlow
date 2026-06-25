nextflow.enable.dsl = 2

process PLOT_GLOBAL_QC {
    tag "Global QC Plotting"
    label 'process_low'


    input:
    path frip_files      
    path cutoff_files    

    output:
    path "FRiPs_*.png"                             , emit: frip_individual_plots
    path "FRiP_comparison_*.png"                   , emit: frip_faceted_plots
    path "FRiP_boxplot.png"                        , emit: frip_boxplot
    path "MACS3_cutoff_peaks_impact_filtered.png"  , emit: cutoff_filtered_plot
    path "Experiment_QC_Summary.csv"               , emit: final_table

    publishDir "${params.outdir}/04_peak_calling/Advanced_QC", mode: 'copy'

    script:
    """
    plot_advanced.qc.R --frip ${frip_files} --cutoffs ${cutoff_files}
    """
}