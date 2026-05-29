nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from '../modules/01_preprocessing/fastqc'
include { FASTQC as FASTQC_TRIM } from '../modules/01_preprocessing/fastqc'
include { TRIM } from '../modules/01_preprocessing/trimming'
include { MULTIQC as MULTIQC_RAW } from '../modules/01_preprocessing/multiqc'
include { MULTIQC as MULTIQC_TRIM } from '../modules/01_preprocessing/multiqc'
include { ALIGNMENT } from '../modules/02_alignment/alignment'
include { PLOT_ALIGNMENT } from '../modules/02_alignment/plotting_alignment'


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

    def trimmed_fastq_ch = Channel.empty() 

    if (params.trim) {
        TRIM(samples_ch)
        trimmed_fastq_ch = TRIM.out
    } else {
            trimmed_fastq_ch = samples_ch.map { sample, f1, f2 ->
                def r1 = file("${params.outdir}/01_qc/trimmed/${sample}/${sample}_trimmed_R1.fastq.gz")
                def r2 = file("${params.outdir}/01_qc/trimmed/${sample}/${sample}_trimmed_R2.fastq.gz")
            
            if ( !r1.exists() || !r2.exists() ) {
                error "Inconsistencia: No existen los archivos recortados para la muestra ${sample} en la carpeta results/. ¿Te has olvidado de activar --trim true?"
            }
            
            return tuple(sample, r1, r2) 
        }
    }

    if (params.fastqc_trim) {
        FASTQC_TRIM(trimmed_fastq_ch, 'fastqc_trimmed')
        MULTIQC_TRIM( Channel.empty().mix(*FASTQC_TRIM.out).collect(), 'trimmed' )    
    }

    if (params.alignment) {
        ALIGNMENT(trimmed_fastq_ch, params.bowtie2_index)
        PLOT_ALIGNMENT( ALIGNMENT.out.summary.collect() )
    }

}