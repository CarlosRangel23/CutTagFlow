nextflow.enable.dsl = 2

process TRIM {

    tag "$sample"

    input:
    tuple val(sample), path(fastq1), path(fastq2)

    output:
    tuple val(sample),
          path("${sample}_trimmed_R1.fastq.gz"),
          path("${sample}_trimmed_R2.fastq.gz")

    publishDir "${params.outdir}/01_qc/trimmed/${sample}",
          mode: 'symlink'


    script:
    """
    fastp \
      -i $fastq1 -I $fastq2 \
      -o ${sample}_trimmed_R1.fastq.gz \
      -O ${sample}_trimmed_R2.fastq.gz \
      -w ${task.cpus} \
      -h ${sample}.fastp.html \
      -j ${sample}.fastp.json \
      -R ${sample} \
      --length_required 30 \
      --trim_poly_g \
      --trim_poly_x \
      --detect_adapter_for_pe

    echo "[TRIM] done for ${sample}"
    """
}
