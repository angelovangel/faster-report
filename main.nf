#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Simple Nextflow wrapper around the aangeloo/faster-report Docker image.
 *
 * Usage:
 *   nextflow run main.nf --fastq_dir /path/to/fastq [other options]
 *
 * flowcell / rundate / basecall / type (platform) are auto-detected from the
 * first fastq file's header (GET_HEADER_DATA). Any of --type / --flowcell /
 * --rundate / --basecall passed explicitly by the user takes precedence over
 * the detected value.
 */

params.fastq     = null          // required: path to folder with fastq files
params.regex     = 'fastq(\\.gz)?$'
params.type      = null          // illumina | ont | pacbio; auto-detected if not set
params.rundate   = null
params.flowcell  = null
params.basecall  = null
params.user      = null
params.save_raw  = false
params.subsample = 1.0
params.outfile   = 'faster-report.html'
params.outdir    = 'output'

if (!params.fastq) {
    error "Please provide a path to a folder with fastq files: --fastq /path/to/fastq"
}

process GET_HEADER_DATA {

    input:
    path fastq

    output:
    path 'header_data.csv'

    script:
    """
    get-header-data.sh ${fastq} "${params.regex}" > header_data.csv
    """
}

process FASTER_REPORT {
    publishDir params.outdir, mode: 'copy'

    input:
    path fastq
    tuple val(platform_detected), val(flowcell_detected), val(rundate_detected), val(basecall_detected)

    output:
    path params.outfile

    script:
    // user-supplied params always win; otherwise fall back to the value
    // detected from the fastq header (GET_HEADER_DATA), unless that is 'NA'
    // (or, for platform, 'unknown')
    def resolve = { userVal, detectedVal, unknownVal ->
        if (userVal) return userVal
        if (detectedVal && detectedVal != unknownVal) return detectedVal
        return null
    }

    def type     = resolve(params.type,     platform_detected, 'unknown') ?: 'ont'
    def rundate  = resolve(params.rundate,  rundate_detected,  'NA')
    def flowcell = resolve(params.flowcell, flowcell_detected, 'NA')
    def basecall = resolve(params.basecall, basecall_detected, 'NA')

    def rundateOpt  = rundate  ? "-d '${rundate}'"  : ''
    def flowcellOpt = flowcell ? "-f '${flowcell}'" : ''
    def basecallOpt = basecall ? "-b '${basecall}'" : ''
    def user        = params.user     ? "-u '${params.user}'" : ''
    def saveraw     = params.save_raw ? "-s TRUE" : ''
    """
    /temp/faster-report.R \\
        -p ${fastq} \\
        -r '${params.regex}' \\
        -t ${type} \\
        ${rundateOpt} \\
        ${flowcellOpt} \\
        ${basecallOpt} \\
        ${user} \\
        ${saveraw} \\
        -x ${params.subsample} \\
        -o ${params.outfile}
    """
}

workflow {
    fastq_ch = Channel.fromPath(params.fastq, type: 'dir', checkIfExists: true)

    header_ch = GET_HEADER_DATA(fastq_ch)
        .splitCsv()
        .map { row -> tuple(row[0], row[1], row[2], row[3]) }
    //header_ch.view()
    FASTER_REPORT(fastq_ch, header_ch)
}
