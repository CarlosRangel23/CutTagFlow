nextflow.enable.dsl = 2

// =========================================================================
// MODULE INCLUDES
// =========================================================================
include { FASTQC as FASTQC_RAW } from '../modules/01_preprocessing/fastqc'
include { FASTQC as FASTQC_TRIM } from '../modules/01_preprocessing/fastqc'
include { TRIM } from '../modules/01_preprocessing/trimming'
include { MULTIQC as MULTIQC_RAW } from '../modules/01_preprocessing/multiqc'
include { MULTIQC as MULTIQC_TRIM } from '../modules/01_preprocessing/multiqc'
include { ALIGNMENT } from '../modules/02_alignment/alignment'
// include { SPIKEIN_ALIGNMENT } from '../modules/02_alignment/spike_in_alignment'
include { PLOT_ALIGNMENT } from '../modules/02_alignment/plotting_alignment'
include { MARK_DUPLICATES } from '../modules/02_alignment/mark_duplicates'
// include { REMOVE_DUPLICATES } from '../modules/02_alignment/remove_duplicates'
include { PLOT_QC_METRICS } from '../modules/02_alignment/plot_alignment_metrics'
// include { PLOT_QC_METRICS as PLOT_QC_FILTERED } from '../modules/02_alignment/plot_alignment_metrics'
// include { PLOT_QC_METRICS as PLOT_QC_FILTERED_DEDUP } from '../modules/02_alignment/plot_alignment_metrics'
// include { FILTER_BAM as FILTER_DUP_BAM } from '../modules/03_filtering/filtering'
// include { FILTER_BAM as FILTER_DEDUP_BAM } from '../modules/03_filtering/filtering'
// include { MACS3 as CALLPEAK_WITH_DUPS } from '../modules/05_peak_calling/macs3'
// include { MACS3 as CALLPEAK_DEDUP } from '../modules/05_peak_calling/macs3'


workflow CUTTAG {

    main:

    // -------------------------------------------------------------------------
    // INITIALIZATION: Parse samplesheet input
    // -------------------------------------------------------------------------
    Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample,
                file(row.fastq1),
                file(row.fastq2),
                row.histone_mark
            )
        }
        .set { samples_ch }

    // Initialize terminal downstream channels to prevent crashes
    def final_filtered_dup_bam_ch   = Channel.empty()
    def final_filtered_dedup_bam_ch = Channel.empty()

    // -------------------------------------------------------------------------
    // PREPROCESSING & QUALITY CONTROL
    // -------------------------------------------------------------------------
    if (params.fastqc) {
        FASTQC_RAW(samples_ch, "fastqc_raw")

        FASTQC_RAW.out
                  .map { zip, html, histone_mark ->
                  return tuple(histone_mark, [zip, html]) }

                  .groupTuple()
                  .map { histone_mark, filesList -> tuple(histone_mark, filesList.flatten()) }
                  .set { raw_fastqc_grouped_ch }

        MULTIQC_RAW( raw_fastqc_grouped_ch, 'raw' )    
    }

    def trimmed_fastq_ch = Channel.empty() 

    if (params.trim) {
        TRIM(samples_ch)
        trimmed_fastq_ch = TRIM.out.map { sample, histone_mark, r1, r2 -> 
            tuple(sample, r1, r2, histone_mark) }
    
    } else {
            trimmed_fastq_ch = samples_ch.map { sample, f1, f2, histone_mark ->            
                def r1 = file("${params.outdir}/01_qc/trimmed/${histone_mark}/${sample}/${sample}_trimmed_R1.fastq.gz")
                def r2 = file("${params.outdir}/01_qc/trimmed/${histone_mark}/${sample}/${sample}_trimmed_R2.fastq.gz")
            
            if ( !r1.exists() || !r2.exists() ) {
                error "No trimmed files for ${sample} (Mark: ${histone_mark}) in results folder"            
            }
            return tuple(sample, r1, r2) 
        }
    }

    if (params.fastqc_trim) {
        FASTQC_TRIM(trimmed_fastq_ch, 'fastqc_trimmed')

        FASTQC_TRIM.out
                   .map { zip, html, histone_mark -> 
                   return tuple(histone_mark, [zip, html]) }

                   .groupTuple()
                   .map { histone_mark, filesList -> tuple(histone_mark, filesList.flatten()) }
                   .set { trim_fastqc_grouped_ch }
    
        MULTIQC_TRIM( trim_fastqc_grouped_ch, 'trimmed' )    
    }

    // -------------------------------------------------------------------------
    // ALIGNMENT & INITIAL QUALITY METRICS
    // -------------------------------------------------------------------------
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

//    if (params.spikein == 'dm6' || params.spikein == 'ecoli') {     
//        SPIKEIN_ALIGNMENT(trimmed_fastq_ch, params.spikein)
//        PLOT_ALIGNMENT( SPIKEIN_ALIGNMENT.out.summary.collect() )        
//    }
//
//    // -------------------------------------------------------------------------
//    // PARALLEL FILTERING BRANCHES (With Duplicates vs Deduplicated)
//    // -------------------------------------------------------------------------
//    if (params.filtering) {
//        def bam_marked_ch = params.alignment ? 
//            MARK_DUPLICATES.out.bam : 
//            samples_ch.map { sample, f1, f2 -> 
//                def file_path = file("${params.outdir}/02_alignment/temp_picard_idxstats/${sample}/${sample}.sorted.dupMarked.bam")
//                if ( !file_path.exists() ) { error "Missing alignment file for ${sample}" }
//                return tuple(sample, file_path)
//            }
//
//        // DUP BRANCH: Preserving Duplicating Traces (Standard practice in CUT&Tag)
//        FILTER_DUP_BAM(bam_marked_ch, 'withDups')
//        PLOT_QC_FILTERED(
//            FILTER_DUP_BAM.out.picard_metrics.collect(),
//            FILTER_DUP_BAM.out.frag_len.collect(),
//            FILTER_DUP_BAM.out.idxstats.collect()
//        )
//        final_filtered_dup_bam_ch = FILTER_DUP_BAM.out.bam
//
//        // DEDUP BRANCH: Purging Duplicating Traces (Control/Comparative branch)
//        REMOVE_DUPLICATES(bam_marked_ch)
//        FILTER_DEDUP_BAM(REMOVE_DUPLICATES.out.bam, 'noDups')
//        PLOT_QC_FILTERED_DEDUP(
//            FILTER_DEDUP_BAM.out.picard_metrics.collect(),
//            FILTER_DEDUP_BAM.out.frag_len.collect(),
//            FILTER_DEDUP_BAM.out.idxstats.collect()
//        )
//        final_filtered_dedup_bam_ch = FILTER_DEDUP_BAM.out.bam
//        
//    } else {
//        // If filtering step is skipped, populate target channels from directory parameters
//        final_filtered_dup_bam_ch = samples_ch.map { sample, f1, f2 ->
//            tuple(sample, file("${params.outdir}/03_filtered/withDups/${sample}.filtered.bam", checkIfExists: params.peaks))
//        }
//        final_filtered_dedup_bam_ch = samples_ch.map { sample, f1, f2 ->
//            tuple(sample, file("${params.outdir}/03_filtered/noDups/${sample}.filtered.dedup.bam", checkIfExists: params.peaks))
//        }
//    }
//
//    // -------------------------------------------------------------------------
//    // DUAL PEAK CALLING LAYER
//    // -------------------------------------------------------------------------
//    if (params.peaks) {
//        CALLPEAK_WITH_DUPS( final_filtered_dup_bam_ch, "noDup" )
//        CALLPEAK_DEDUP( final_filtered_dedup_bam_ch, "deDup" )
//    }
}