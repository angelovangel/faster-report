#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Simple Nextflow wrapper around the aangeloo/faster-report Docker image.
 *
 * Usage:
 *   nextflow run main.nf --fastq_dir /path/to/fastq [other options]
 *
 * flowcell / rundate / basecall are auto-detected from the first fastq
 * file's header (GET_HEADER_DATA). Any of --flowcell / --rundate / --basecall
 * passed explicitly by the user takes precedence over the detected value.
 */

params.fastq_dir = null          // required: path to folder with fastq files
params.regex     = 'fast(q|q.gz)$'
params.type      = 'ont'         // illumina | ont | pacbio
params.rundate   = null
params.flowcell  = null
params.basecall  = null
params.user      = null
params.save_raw  = false
params.subsample = 1.0
params.outfile   = 'faster-report.html'
params.outdir    = 'output'

if (!params.fastq_dir) {
    error "Please provide a path to a folder with fastq files: --fastq_dir /path/to/fastq"
}

process GET_HEADER_DATA {
    tag "header info for ${fastq_dir}"

    input:
    path fastq_dir

    output:
    path 'header_data.csv'

    script:
    """
    get-header-data.sh ${fastq_dir} '${params.regex}' > header_data.csv
    """
}

process FASTER_REPORT {
    publishDir params.outdir, mode: 'copy'

    input:
    path fastq_dir
    tuple val(flowcell_detected), val(rundate_detected), val(basecall_detected)

    output:
    path params.outfile

    script:
    // user-supplied params always win; otherwise fall back to the value
    // detected from the fastq header (GET_HEADER_DATA), unless that is 'NA'
    def resolve = { userVal, detectedVal ->
        if (userVal) return userVal
        if (detectedVal && detectedVal != 'NA') return detectedVal
        return null
    }

    def rundate  = resolve(params.rundate,  rundate_detected)
    def flowcell = resolve(params.flowcell, flowcell_detected)
    def basecall = resolve(params.basecall, basecall_detected)

    def rundateOpt  = rundate  ? "-d '${rundate}'"  : ''
    def flowcellOpt = flowcell ? "-f '${flowcell}'" : ''
    def basecallOpt = basecall ? "-b '${basecall}'" : ''
    def user        = params.user     ? "-u '${params.user}'" : ''
    def saveraw     = params.save_raw ? "-s TRUE" : ''
    """
    /temp/faster-report.R \\
        -p ${fastq_dir} \\
        -r '${params.regex}' \\
        -t ${params.type} \\
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
    fastq_ch = Channel.fromPath(params.fastq_dir, type: 'dir', checkIfExists: true)

    header_ch = GET_HEADER_DATA(fastq_ch)
        .splitCsv()
        .map { row -> tuple(row[0], row[1], row[2]) }

    FASTER_REPORT(fastq_ch, header_ch)
}
