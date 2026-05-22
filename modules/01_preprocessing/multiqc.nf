nextflow.enable.dsl = 2

process MULTIQC {

    label 'qc'

    input:
    path qc_dir

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    publishDir "${params.outdir}/01_qc/multiqc",
               mode: 'copy'

    script:
    """
    multiqc $qc_dir
    """
}
