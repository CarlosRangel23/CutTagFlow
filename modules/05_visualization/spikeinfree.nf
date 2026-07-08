nextflow.enable.dsl = 2

process SPIKE_IN_FREE {
    label 'process_medium'

    publishDir "${params.outdir}/05_visualization/spikein_factors", mode: 'copy'

    input:
    path bams 
    path bais 
    path meta_file
    
    output:
    path "cuttag_factors_SF.txt"           , emit: scaling_factors
    path "cuttag_factors_distribution.pdf" , emit: test_distribution, optional: true
    path "cuttag_factors_boxplot.pdf"      , emit: boxplot, optional: true
    path "cuttag_factors_rawCounts.txt"    , emit: rawCounts
    path "cuttag_factors_parsedMatrix.txt" , emit: cuttag_factors_parsedMatrix

    script:
    """
    Rscript -e "
    library(ChIPseqSpikeInFree)

    bams <- list.files(pattern = '\\\\.bam\$')

    ChIPseqSpikeInFree(
        bamFiles = bams, 
        chromFile = 'hg38', 
        metaFile = '$meta_file', 
        prefix = 'cuttag_factors'
    )
    "
    """
}