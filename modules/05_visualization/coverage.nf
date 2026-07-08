nextflow.enable.dsl = 2

process COVERAGE {
    tag "${sample}"
    label 'process_medium'

    publishDir "${params.outdir}/05_visualization/coverage/${histone_mark}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai), val(histone_mark), 

    output:
    tuple val(sample), path("*.bw"), emit: bigwig

    script:
    """
    Rscript -e "
    library(ChIPseqSpikeInFree)
    write.lines('${bam}', 'bam_list.txt')
    ChIPseqSpikeInFree(bamIndexFile = 'bam_list.txt', chromFile = 'hg38')
    "

    SCALE_FACTOR=\$(awk 'NR==2 {print \$2}' *export.txt || awk 'NR==2 {print \$2}' *_ChIPseqSpikeInFree.txt)

    bamCoverage \\
        -b ${bam} \\
        -o ${sample}.${histone_mark}.normalized.bw \\
        --scaleFactor \$SCALE_FACTOR \\
        -p max \\
        --binSize 10
    """
}