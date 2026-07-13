#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Simple Nextflow wrapper around the aangeloo/faster-report Docker image.
 *
 * Usage:
 *   nextflow run main.nf --fastq /path/to/fastq [other options]
 *
 * flowcell / rundate / basecall / type (platform) are auto-detected from the
 * first fastq file's header (GET_HEADER_DATA). Any of --type / --flowcell /
 * --rundate / --basecall passed explicitly by the user takes precedence over
 * the detected value.
 */

params.reads     = null          // required: path to folder with fastq files
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

if (!params.reads) {
    error "Please provide a path to a folder with fastq/bam files: --reads /path/to/fastq"
}

process CONVERT_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    
    input:
        path reads

    output:
        path("*.fastq")

    script:
    """
    samtools fastq -@ ${task.cpus} -T '*' ${reads} > ${reads.simpleName}.fastq
    """

}

process GET_HEADER_DATA {

    input:
    path 'temp/*' // collected there and passed to script

    output:
    path 'header_data.csv'

    script:
    """
    get-header-data.sh temp "${params.regex}" > header_data.csv
    """
}

process FASTER_REPORT {
    publishDir params.outdir, mode: 'copy'

    input:
    path 'temp/*'
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
        -p temp \\
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

pattern = "*.{bam,fasta,fastq,fastq.gz,fq,fq.gz}"
ch_files = Channel.fromPath(params.reads + "/" + pattern, type: 'file', checkIfExists: true) 

workflow {
    // Branch reads into BAM and others for centralized conversion
    //ch_reads = Channel.fromPath(params.reads, type: 'dir', checkIfExists: true)
    ch_files
        .branch {
            bam: it.name.endsWith('.bam')
            other: true
        }
        .set { ch_reads_split }

    ch_fastq = CONVERT_READS(ch_reads_split.bam)
        .mix(ch_reads_split.other)
        .collect()

    //ch_fastq.view()
    //fastq_ch = Channel.fromPath(params.reads, type: 'dir', checkIfExists: true)

    header_ch = GET_HEADER_DATA(ch_fastq)
        .splitCsv()
        .map { row -> tuple(row[0], row[1], row[2], row[3]) }
    //header_ch.view()
    FASTER_REPORT(ch_fastq, header_ch)
}
