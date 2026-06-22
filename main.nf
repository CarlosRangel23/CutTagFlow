#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { CUTTAG } from './workflows/cuttag.nf'


workflow {
    CUTTAG()
}