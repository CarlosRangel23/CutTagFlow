nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from '../modules/01_preprocessing/fastqc'
include { FASTQC as FASTQC_TRIM } from '../modules/01_preprocessing/fastqc'
include { TRIM } from '../modules/01_preprocessing/trimming'
include { MULTIQC as MULTIQC_RAW } from '../modules/01_preprocessing/multiqc'
include { MULTIQC as MULTIQC_TRIM } from '../modules/01_preprocessing/multiqc'
include { ALIGNMENT } from '../modules/02_alignment/alignment'


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
        FASTQC_RAW(samples_ch,"fastqc_raw")
        MULTIQC_RAW( Channel.empty().mix(*FASTQC_RAW.out).collect(), 'raw' )    
    }

    if (params.trim) {
        TRIM(samples_ch)
    }
    
    def trimmed_fastq_ch = params.trim ? 
        TRIM.out : 
        samples_ch.map { sample, f1, f2 ->
            tuple(
                sample,
                file("${params.outdir}/01_qc/trimmed/${sample}/${sample}_trimmed_R1.fastq.gz", checkIfExists: true),
                file("${params.outdir}/01_qc/trimmed/${sample}/${sample}_trimmed_R2.fastq.gz", checkIfExists: true)
            )
        }

    if (params.fastqc_trim) {
        FASTQC_TRIM(trimmed_fastq_ch, 'fastqc_trimmed')
        MULTIQC_TRIM( Channel.empty().mix(*FASTQC_TRIM.out).collect(), 'trimmed' )    
    }

    if (params.alignment) {
        ALIGNMENT(trimmed_fastq_ch, params.bowtie2_index)
    }

}