nextflow.enable.dsl = 2

process SPIKE_IN_FREE {
    tag "All_Samples"
    label 'process_medium'

    publishDir "${params.outdir}/05_visualization/spikein_factors", mode: 'copy'

    input:
    path bams 
    path bais 
    path meta_file
    val chromFile
    
    output:
    path "cuttag_factors_SF.txt"           , emit: scaling_factors
    path "cuttag_factors_distribution.pdf" , emit: test_distribution, optional: true
    path "cuttag_factors_boxplot.pdf"      , emit: boxplot, optional: true
    path "cuttag_factors_rawCounts.txt"    , emit: rawCounts
    path "cuttag_factors_parsedMatrix.txt" , emit: parsedMatrix

    script:
    """
    #!/usr/bin/env Rscript
    library(ChIPseqSpikeInFree)
    
    meta <- read.table("${meta_file}", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    bams_vector <- as.character(meta\$ID)

    ChIPseqSpikeInFree(
        bamFiles = bams_vector, 
        chromFile = "${chromFile}", 
        metaFile = "${meta_file}", 
        prefix = "cuttag_factors"
    )
    """
}