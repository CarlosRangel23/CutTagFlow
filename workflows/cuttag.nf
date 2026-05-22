nextflow.enable.dsl = 2

include { FASTQC } from '../modules/01_preprocessing/fastqc'
include { TRIM } from '../modules/01_preprocessing/trimming'

workflow CUTTAG {

    main:

    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample,
                file(row.fastq1),
                file(row.fastq2)
            )
        }
        .set { samples_ch }


    if (params.fastqc) {
        FASTQC(samples_ch,"fastqc_raw")
    }

    if (params.trim) {
        TRIM(samples_ch)
            reads_for_alignment = TRIM.out.reads
    }

    if (params.fastqc_trim) {
        samples_ch
            .map { sample, fastq1, fastq2 ->
                tuple(
                    sample,
                    file("${params.outdir}/01_qc/trimmed/${sample}/${sample}_trimmed_R1.fastq.gz"),
                    file("${params.outdir}/01_qc/trimmed/${sample}/${sample}_trimmed_R2.fastq.gz")
                )
            }
            .set { trimmed_fastq_ch }
    
        FASTQC(trimmed_fastq_ch, 'fastqc_trimmed')
    }

    if (params.alignment) {
        ALIGNMENT(samples_ch)
    }

    if (params.peaks) {
        PEAKS(samples_ch)
    }

}