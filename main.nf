#! /usr/bin/env nextflow

nextflow.enable.dsl = 2




include { flexAMR        } from './workflows/flexAMR'

workflow {
	flexAMR ()
	
}


	



