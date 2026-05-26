nextflow.enable.dsl = 2

process FASTQC {

    tag "$sample"
    label 'qc'

    input:
    tuple val(sample), path(fastq1), path(fastq2)
    val qc_label

    output:
    path "*_fastqc.zip"
    path "*_fastqc.html"

    publishDir "${params.outdir}/01_qc/${qc_label}/${sample}",
          mode: 'copy'

    script:
    """
    fastqc \
      --t ${task.cpus} \
      --o . \
      $fastq1 $fastq2

    echo "[FASTQC] ${qc_label} QC done for ${sample}"
    """
}
