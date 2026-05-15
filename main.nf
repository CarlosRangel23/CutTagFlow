#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUTTAG } from './workflows/cuttag.nf'

workflow {
    if (params.all) {
    params.fastqc    = true
    params.trim      = true
    params.fastqc_trim = true
    params.alignment = true
    params.peaks     = true
}

    CUTTAG()
}
