nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW } from '../modules/01_preprocessing/fastqc'
include { FASTQC as FASTQC_TRIM } from '../modules/01_preprocessing/fastqc'
include { TRIM } from '../modules/01_preprocessing/trimming'
include { MULTIQC as MULTIQC_RAW } from '../modules/01_preprocessing/multiqc'
include { MULTIQC as MULTIQC_TRIM } from '../modules/01_preprocessing/multiqc'
include { ALIGNMENT } from '../modules/02_alignment/alignment'
include { SPIKEIN_ALIGNMENT } from '../modules/02_alignment/spike_in_alignment'
include { PLOT_ALIGNMENT } from '../modules/02_alignment/plotting_alignment'
include { MARK_DUPLICATES } from '../modules/02_filtering_conversion/mark_duplicates'
include { PLOT_QC_METRICS } from '../modules/02_filtering_conversion/plot_alignment_metrics'


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
                error "No trimmed files for ${sample} in results folder"
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
        MARK_DUPLICATES( ALIGNMENT.out.bam )
        PLOT_QC_METRICS( 
            MARK_DUPLICATES.out.picard_metrics.collect(), 
            MARK_DUPLICATES.out.frag_len.collect(), 
            MARK_DUPLICATES.out.idxstats.collect() 
        )    
    }

    if (params.spikein == 'dm6' || params.spikein == 'ecoli') {     
        SPIKEIN_ALIGNMENT(trimmed_fastq_ch, params.spikein)
        PLOT_ALIGNMENT( SPIKEIN_ALIGNMENT.out.summary.collect() )        
    }

}