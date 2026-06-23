nextflow.enable.dsl = 2

process MACS3 {
    tag "MACS3 Callpeak (${sample} - ${label})"
    label 'process_medium'

    input:
    tuple val(sample), path(bam), path(bai), val(histone_mark)
    val label
    val genome_size

    output:
    path "${sample}.${label}_peaks.*"                                     , emit: peaks
    path "${sample}.${label}_summits.bed"                                 , emit: summits, optional: true
    tuple val(sample), val(histone_mark), path("*.{narrowPeak,broadPeak}"), emit: peak_file
    path "*_cutoff_analysis.txt"                                          , emit: summary

    publishDir "${params.outdir}/04_peak_calling/${label}/${histone_mark}/${sample}", mode: 'copy'

    script:
    def macs3_args = ""
    if ( histone_mark in ['H3K27me3', 'H3K9me3', 'H3K36me3', 'H3K4me1', 'H3K4me2'] ) {
        macs3_args = "--broad --broad-cutoff 1e-5 --nolambda"
    } else {
        macs3_args = "-q 1e-5 --nolambda"
    }

    // References for the parameters:
    // https://macs3-project.github.io/MACS/docs/callpeak.html
    // Kaya-Okur HS et al., Nat Commun, 2019 (Original CUT&Tag protocol description)
    // Cheng S et al., Genomics Proteomics Bioinformatics, 2024 (Tool performance benchmarking)
    // Abbasova L et al., Nat Commun, 2025 (MACS parameter optimization for CUT&Tag against ENCODE)

    """
    macs3 callpeak \\
        -t ${bam} \\
        -f BAMPE \\
        -g ${genome_size} \\
        -n ${sample}.${label} \\
        ${macs3_args} \\
        --nomodel \\
        --keep-dup all \\
        --cutoff-analysis
    """
}