nextflow.enable.dsl = 2

process TSSE {
    tag "${sample} - ${label}"
    label 'process_medium' // Requiere más RAM para cargar el BAM en R y desplazar los reads

    input:
    tuple val(sample), val(histone_mark), path(bam), path(bai), path(peak)
    val label

    output:
    tuple val(sample), val(histone_mark), val(label), path("${sample}.${label}.tsse.csv"), emit: tsse_csv

    publishDir "${params.outdir}/04_peak_calling/${label}/${histone_mark}/TSSE", mode: 'copy'
    
    script:
    """
    #!/usr/bin/env Rscript

    # 1. Load libraries
    suppressPackageStartupMessages({
      library(ATACseqQC)
      library(Rsamtools)
      library(GenomicAlignments)
      library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    })

    # 2. Define gene transcripts for TSS profiling
    txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

    # 3. Read BAM file with paired-end configuration
    gal <- readBamFile(
      "${bam}",
      asMates = TRUE,
      bigFile = TRUE
    )

    # 4. Shift Tn5 insertions and write temporary shifted BAM
    shiftedBamfile <- "${sample}_shifted.bam"
    Gal1 <- shiftGAlignmentsList(
      gal,
      outbam = shiftedBamfile
    )

    # 5. Index the newly created shifted BAM file
    indexBam(shiftedBamfile)

    # 6. Calculate TSSE score using the shifted alignments
    tsse_result <- TSSeScore(Gal1, txs)
    tsse_value  <- tsse_result$TSSeScore

    # 7. Consolidate into a standardized CSV for Nextflow matching
    tsse_data <- data.frame(
      sample       = "${sample}",
      histone_mark = "${histone_mark}",
      label        = "${label}",
      tsse_score   = tsse_value,
      stringsAsFactors = FALSE
    )

    write.csv(tsse_data, file = "${sample}.${label}.tsse.csv", row.names = FALSE)

    # 8. Clean up heavy intermediate BAM files from the work directory
    file.remove(shiftedBamfile)
    file.remove(paste0(shiftedBamfile, ".bai"))
    """
}