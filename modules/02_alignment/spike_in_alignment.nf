nextflow.enable.dsl = 2

process SPIKEIN_ALIGNMENT {

    tag "$sample"
    label 'process_high' 

    publishDir "${params.outdir}/02_alignment/spike_in/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(fastq1), path(fastq2)
    val spikein_organism 
    
    output:
    tuple val(sample), path("${sample}.sorted.spikein.bam"), emit: bam
    path "${sample}.sorted.spikein.bam.bai"                 , emit: bai
    path "${sample}.spikein.bowtie2.txt"                    , emit: summary

    script:
    def cores = task.cpus
    def ref = "data/reference/${spikein_organism}/${spikein_organism}"

    """
    bowtie2 \
      --end-to-end \
      --very-sensitive \
      --no-mixed \
      --no-discordant \
      --phred33 \
      -I 10 \
      -X 700 \
      -p ${cores} \
      -x ${ref} \
      -1 ${fastq1} \
      -2 ${fastq2} \
      --rg-id ${sample} \
      --rg "SM:${sample}" \
      --rg "LB:lib_${sample}" \
      --rg "PL:ILLUMINA" \
    | samtools sort -@ ${cores} -o ${sample}.sorted.spikein.bam -

    samtools index ${sample}.sorted.spikein.bam

    cp .command.log ${sample}.spikein.bowtie2.txt || true
    
    echo "[SPIKE-IN] Aligned spike-in bam generated for ${sample} against ${spikein_organism}"
    """
}