nextflow.enable.dsl = 2

process MULTIQC {

    label 'qc'

    input:
    path qc_files
    val qc_label

    output:
    path "${qc_label}_multiqc_report.html"
    path "${qc_label}_multiqc_report_data"

    publishDir "${params.outdir}/01_qc/multiqc/${qc_label}", mode: 'copy'

    script:
    """
    multiqc . -n ${qc_label}_multiqc_report.html
    """
}