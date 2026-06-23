nextflow.enable.dsl = 2

process FRIP_SCORE {
    tag "${sample} - ${label}"
    label 'process_low'

    input:
    tuple val(sample), val(histone_mark), path(bam), path(bai), path(peak)
    val label
    path tss_bed
    path chr_sizes

    output:
    tuple val(sample), val(histone_mark), val(label), path("${sample}.${label}.frip.csv"), emit: metrics_csv
    
    publishDir "${params.outdir}/04_peak_calling/${label}/${histone_mark}/FRIP", mode: 'copy'


 script:
    """
    # 1. Expand TSS regions ±2 kb using bedtools slop
    bedtools slop \\
        -i ${tss_bed} \\
        -g ${chr_sizes} \\
        -b 2000 \\
        > tss_2kb.bed

    # 2. Count overlapping reads directly into Bash variables
    reads_tss=\$(bedtools intersect -a ${bam} -b tss_2kb.bed -u -bed | wc -l)
    reads_peaks=\$(bedtools intersect -a ${bam} -b ${peak} -u -bed | wc -l)

    # 3. Total mapped reads
    total_reads=\$(samtools view -c ${bam})

    # 4. Perform arithmetic directly in Bash using bc
    if [ "\$total_reads" -gt 0 ]; then
        frip_tss=\$(echo "scale=6; \$reads_tss / \$total_reads" | bc -l)
        frip_peaks=\$(echo "scale=6; \$reads_peaks / \$total_reads" | bc -l)
    else
        frip_tss=0
        frip_peaks=0
    fi

    # 5. Build the CSV natively by echoing the header and rows
    echo "sample,histone_mark,label,total_reads,frip_tss_2kb,frip_peaks" > ${sample}.${label}.frip.csv
    echo "${sample},${histone_mark},${label},\$total_reads,\$frip_tss,\$frip_peaks" >> ${sample}.${label}.frip.csv

    # 6. Remove remaining intermediate files
    rm tss_2kb.bed
    """
}