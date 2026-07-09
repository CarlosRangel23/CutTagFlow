nextflow.enable.dsl = 2

process COVERAGE {
    tag "${sample}"
    label 'process_medium'

    publishDir "${params.outdir}/05_visualization/coverage/${histone_mark}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(bam), path(bai), val(histone_mark), path (sf_file)

    output:
    tuple val(sample), path("*.bw"), emit: bigwig

    script:
    """
    SCALE_FACTOR=\$(grep "${bam.name}" ${sf_file} | awk '{print \$7}')
    
    bamCoverage \\
        -b ${bam} \\
        -o ${sample}.bw \\
        -of "bigwig" \\
        -e \\
        -p ${task.cpus} \\
        --scaleFactor \$SCALE_FACTOR \\
        --binSize 10
    """
}