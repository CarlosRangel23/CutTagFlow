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
include { PLOT_ALIGNMENT } from '../modules/02_alignment/plotting_alignment'
include { MARK_DUPLICATES } from '../modules/02_alignment/mark_duplicates'
include { PLOT_QC_METRICS } from '../modules/02_alignment/plot_alignment_metrics'
include { PLOT_QC_FILTERED as PLOT_QC_FILTERED } from '../modules/03_filtering/plot_filtered_metrics'
include { PLOT_QC_FILTERED as PLOT_QC_FILTERED_DEDUP } from '../modules/03_filtering/plot_filtered_metrics'
include { FILTER_BAM as FILTER_DUP_BAM } from '../modules/03_filtering/filtering'
include { FILTER_BAM as FILTER_DEDUP_BAM } from '../modules/03_filtering/filtering'
include { MACS3 as CALLPEAK_DEDUP } from '../modules/04_peak_calling/macs3'
include { MACS3 as CALLPEAK_WITH_DUPS } from '../modules/04_peak_calling/macs3'
include { FRIP_SCORE as FRIP_SCORE_DUP } from '../modules/04_peak_calling/frip_score'
include { FRIP_SCORE as FRIP_SCORE_DEDUP } from '../modules/04_peak_calling/frip_score'
include { PLOT_GLOBAL_QC } from '../modules/04_peak_calling/plot_global_qc'
include { CONSENSUS } from '../modules/04_peak_calling/consensus'
include { DIFFBIND } from '../modules/04_peak_calling/diffbind'
include { SPIKE_IN_FREE } from '../modules/05_visualization/spikeinfree'
include { COVERAGE as DEEPTOOLS_COVERAGE_DEDUPS } from '../modules/05_visualization/coverage'
include { GATK_HAPLOTYPE_CALLER } from '../modules/06_variant_calling/haplotype_caller'
include { GATK_JOINT_GENOTYPING } from '../modules/06_variant_calling/joint_genotyping'


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
            return tuple(sample, r1, r2, histone_mark) 
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

    // -------------------------------------------------------------------------
    // PARALLEL FILTERING BRANCHES (With Duplicates vs Deduplicated)
    // -------------------------------------------------------------------------
    if (params.filtering) {
        def bam_marked_ch = params.alignment ? 
            MARK_DUPLICATES.out.bam.join(MARK_DUPLICATES.out.bai)
            .map { sample, bam, histone_mark, bai ->
                return tuple(sample, bam, bai, histone_mark) 
            } : 
            samples_ch.map { sample, f1, f2, histone_mark -> 
                def bam_path = file("${params.outdir}/02_alignment/temp_picard_idxstats/${sample}/${sample}.sorted.dupMarked.bam")
                def bai_path = file("${params.outdir}/02_alignment/temp_picard_idxstats/${sample}/${sample}.sorted.dupMarked.bam.bai")
                if ( !bam_path.exists() ) { error "Missing alignment file for ${sample}" }
                return tuple(sample, bam_path, bai_path, histone_mark)
            }

        // DUP BRANCH: Preserving Duplicating Traces (Standard practice in CUT&Tag)
        // FROM HENIKOFF LAB: CUT&Tag integrates adapters into DNA in the vicinity of the antibody-tethered pA-Tn5, and the exact sites of integration are affected 
        // by the accessibility of surrounding DNA. For this reason fragments that share exact starting and ending positions are expected to be common, and such ‘duplicates’
        // may not be due to duplication during PCR. In practice, we have found that the apparent duplication rate is low for high quality CUT&Tag datasets, and even the
        // apparent ‘duplicate’ fragments are likely to be true fragments. Thus, we DO NOT recommend removing the duplicates. In experiments with very small amounts of 
        // material or where PCR duplication is suspected, duplicates can be removed.

        FILTER_DUP_BAM(bam_marked_ch, 'withDups', params.blacklist_bed)
        PLOT_QC_FILTERED(
            FILTER_DUP_BAM.out.frag_len.collect(),
            FILTER_DUP_BAM.out.idxstats.collect(),
            'withDups'
        )

        final_filtered_dup_bam_ch = FILTER_DUP_BAM.out.bam
                                                  .join(FILTER_DUP_BAM.out.bai)
                                                  .map { sample, bam, label, histone_mark, bai  ->
                                                    return tuple(sample, bam, bai, histone_mark) 
                                                  }

        // DEDUP BRANCH: Purging Duplicating Traces (Control/Comparative branch)
        FILTER_DEDUP_BAM(bam_marked_ch, 'noDups', params.blacklist_bed)
        PLOT_QC_FILTERED_DEDUP(
            FILTER_DEDUP_BAM.out.frag_len.collect(),
            FILTER_DEDUP_BAM.out.idxstats.collect(),
            'noDups'
        )

        final_filtered_dedup_bam_ch = FILTER_DEDUP_BAM.out.bam
                                                  .join(FILTER_DEDUP_BAM.out.bai)
                                                  .map { sample, bam, label, histone_mark, bai  ->
                                                    return tuple(sample, bam, bai, histone_mark) 
                                                  }
        
    } else {
        // If filtering step is skipped, populate target channels from directory parameters
        final_filtered_dup_bam_ch = samples_ch.map { sample, f1, f2, histone_mark ->
            def bam_path = file("${params.outdir}/03_filtered/${histone_mark}/withDups/${sample}/${sample}.withDups.filtered.bam", checkIfExists: params.peaks)
            def bai_path = file("${params.outdir}/03_filtered/${histone_mark}/withDups/${sample}/${sample}.withDups.filtered.bam.bai", checkIfExists: params.peaks)
            return tuple(sample, bam_path, bai_path, histone_mark)        
                }

        final_filtered_dedup_bam_ch = samples_ch.map { sample, f1, f2, histone_mark ->
            def bam_path = file("${params.outdir}/03_filtered/${histone_mark}/noDups/${sample}/${sample}.noDups.filtered.bam", checkIfExists: params.peaks)
            def bai_path = file("${params.outdir}/03_filtered/${histone_mark}/noDups/${sample}/${sample}.noDups.filtered.bam.bai", checkIfExists: params.peaks)
            return tuple(sample, bam_path, bai_path, histone_mark)      
                }
    }

    // -------------------------------------------------------------------------
    // DUAL PEAK CALLING LAYER
    // -------------------------------------------------------------------------
    if (params.peaks) {
        CALLPEAK_WITH_DUPS(final_filtered_dup_bam_ch, "withDups", params.genome_size)
        CALLPEAK_DEDUP(final_filtered_dedup_bam_ch, "noDups", params.genome_size)
    

        // -------------------------------------------------------------------------
        // ADVANCED QC
        // -------------------------------------------------------------------------

        dup_qc_input_ch = final_filtered_dup_bam_ch.join(CALLPEAK_WITH_DUPS.out.peak_file)
                                                   .map { sample, bam, bai, histone_mark, histone_mark_rep, peak ->
                                                   return tuple(sample, histone_mark, bam, bai, peak) }
    
        dedup_qc_input_ch = final_filtered_dedup_bam_ch.join(CALLPEAK_DEDUP.out.peak_file)
                                                       .map { sample, bam, bai, histone_mark, histone_mark_rep, peak ->
                                                       return tuple(sample, histone_mark, bam, bai, peak) }

        FRIP_SCORE_DUP(dup_qc_input_ch,"withDups", params.tss_bed, params.chr_sizes)
        FRIP_SCORE_DEDUP(dedup_qc_input_ch,"noDups", params.tss_bed, params.chr_sizes)

        all_qc_files_ch = FRIP_SCORE_DUP.out.metrics_csv
                                .mix(FRIP_SCORE_DEDUP.out.metrics_csv)
                                .flatMap { sample, mark, label, frip_csv -> [frip_csv] }
                                .collect()

        all_cutoffs_ch = CALLPEAK_WITH_DUPS.out.summary
                                           .mix(CALLPEAK_DEDUP.out.summary)
                                           .collect()

        PLOT_GLOBAL_QC(all_qc_files_ch, all_cutoffs_ch) 

        //-------------------------------------------------------------------------
        // CONSENSUS PEAKS
        // -------------------------------------------------------------------------
        
        with_dups_peaks = CALLPEAK_WITH_DUPS.out.peak_file
                                            .map { sample, histone_mark, peak -> 
                                            tuple( tuple("withDups", histone_mark), peak ) }

        no_dups_peaks = CALLPEAK_DEDUP.out.peak_file
                                      .map { sample, histone_mark, peak ->
                                      tuple( tuple("noDups", histone_mark), peak ) }

        with_dups_bams = final_filtered_dup_bam_ch
                                            .map { sample, bam, bai, histone_mark -> tuple( tuple("withDups", histone_mark), bam ) }

        no_dups_bams = final_filtered_dedup_bam_ch
                                            .map { sample, bam, bai, histone_mark -> tuple( tuple("noDups", histone_mark), bam ) }        

        all_peaks_ch = with_dups_peaks.mix(no_dups_peaks).groupTuple()
        all_bams_ch  = with_dups_bams.mix(no_dups_bams).groupTuple()
        consensus_input_ch = all_peaks_ch.join(all_bams_ch)

        CONSENSUS(consensus_input_ch)
        
    }

    if (params.diffbind){
        if (params.peaks){
            diffbind_input_ch = consensus_input_ch.map { meta, peaks, bams ->
                                                tuple(meta, peaks, bams, params.min_overlap) }
            DIFFBIND(diffbind_input_ch)
        }
        else {
            with_dups_fallback = samples_ch.map { sample, f1, f2, histone_mark ->
                def peak_pattern = "${params.outdir}/04_peak_calling/withDups/${histone_mark}/${sample}/${sample}.withDups_peaks.{narrowPeak,broadPeak}"
                def peak_files   = files(peak_pattern)
                def peak_path    = peak_files ? peak_files[0] : error("No peak file found for ${sample} (withDups) in results directory.")
                def bam_path     = file("${params.outdir}/03_filtered/${histone_mark}/withDups/${sample}/${sample}.withDups.filtered.bam", checkIfExists: true)

                return tuple( tuple("withDups", histone_mark), peak_path, bam_path )
            }

            no_dups_fallback = samples_ch.map { sample, f1, f2, histone_mark ->
                def peak_pattern = "${params.outdir}/04_peak_calling/noDups/${histone_mark}/${sample}/${sample}.noDups_peaks.{narrowPeak,broadPeak}"
                def peak_files   = files(peak_pattern)
                def peak_path    = peak_files ? peak_files[0] : error("No peak file found for ${sample} (noDups) in results directory.")
                def bam_path     = file("${params.outdir}/03_filtered/${histone_mark}/noDups/${sample}/${sample}.noDups.filtered.bam", checkIfExists: true)
                
                return tuple( tuple("noDups", histone_mark), peak_path, bam_path )
            }

            all_fallback_peaks_ch = with_dups_fallback.mix(no_dups_fallback)
                .map { meta, peak, bam -> tuple(meta, peak) }
                .groupTuple()

            all_fallback_bams_ch = with_dups_fallback.mix(no_dups_fallback)
                .map { meta, peak, bam -> tuple(meta, bam) }
                .groupTuple()

            diffbind_input_ch = all_fallback_peaks_ch.join(all_fallback_bams_ch)
                .map { meta, peaks, bams -> tuple(meta, peaks, bams, params.min_overlap) }

            DIFFBIND(diffbind_input_ch)
        }    
    }

    // ---------------------------------------------------------------------
    // QUANTITATIVE VISUALIZATION GENERATION (DeepTools)
    // ---------------------------------------------------------------------
    if (params.coverage) {
        final_filtered_dedup_bam_ch
            .multiMap { sample, bam, bai, histone_mark ->
                bams: bam
                bais: bai
            }
            .set { gathered_files_ch }

        SPIKE_IN_FREE( 
            gathered_files_ch.bams.collect(), 
            gathered_files_ch.bais.collect(), 
            file(params.meta_spike), 
            params.chromFile 
        )


        deeptools_input_ch = final_filtered_dedup_bam_ch.combine(SPIKE_IN_FREE.out.scaling_factors)
        
        DEEPTOOLS_COVERAGE_DEDUPS( deeptools_input_ch )
     }

    // ---------------------------------------------------------------------
    // VARIANT CALLING 
    // ---------------------------------------------------------------------
    if (params.variant) {
        GATK_HAPLOTYPE_CALLER( 
            final_filtered_dedup_bam_ch, 
            file(params.genome_fasta), 
            file(params.genome_fai), 
            file(params.genome_dict)
        )

        joint_input_ch = GATK_HAPLOTYPE_CALLER.out.gvcf.groupTuple()

        GATK_JOINT_GENOTYPING( 
            joint_input_ch, 
            file(params.genome_fasta), 
            file(params.genome_fai), 
            file(params.genome_dict)
        )
    }

}


