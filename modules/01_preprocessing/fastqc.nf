nextflow.enable.dsl = 2

process FASTQC {

    tag "${histone_mark}/${sample}"
    label 'qc'

    input:
    tuple val(sample), path(fastq1), path(fastq2), val(histone_mark)
    val qc_label

    output:
    tuple path("*_fastqc.zip"), path("*_fastqc.html"), val(histone_mark)

    publishDir "${params.outdir}/01_qc/${qc_label}/${histone_mark}/${sample}",
          mode: 'copy'

    script:
    """
    fastqc \
      --t ${task.cpus} \
      --o . \
      $fastq1 $fastq2

    echo "[FASTQC] ${qc_label} QC done for ${sample} (${histone_mark})"
    """
}
