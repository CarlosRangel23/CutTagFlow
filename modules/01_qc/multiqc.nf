nextflow.enable.dsl = 2

process MULTIQC {

    tag "${qc_label}/${histone_mark}"
    label 'qc'

    input:
    tuple val(histone_mark), path(qc_files)
    val qc_label

    output:
    path "multiqc_${qc_label}_${histone_mark}.html"
    path "multiqc_${qc_label}_${histone_mark}_data"

    publishDir "${params.outdir}/01_qc/multiqc/${qc_label}/${histone_mark}", mode: 'copy'
    
    script:
    """
    multiqc . \
      -n "multiqc_${qc_label}_${histone_mark}" \
      --filename "multiqc_${qc_label}_${histone_mark}.html"
    """
}
