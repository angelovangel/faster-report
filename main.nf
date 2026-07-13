#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Simple Nextflow wrapper around the aangeloo/faster-report Docker image.
 *
 * Usage:
 *   nextflow run main.nf --reads /path/to/fastq (or bam) [other options]
 *
 * flowcell / rundate / basecall / type (platform) are auto-detected from the
 * first fastq file's header (GET_HEADER_DATA). Any of --type / --flowcell /
 * --rundate / --basecall passed explicitly by the user takes precedence over
 * the detected value.
 */

if (params.help) {
    log.info """
    ================================================================================
     FASTER-REPORT PIPELINE
     https://github.com/angelovangel/faster-report
    ================================================================================
     Generates interactive HTML quality control reports for FASTQ/BAM files.

     Usage:
       nextflow run angelovangel/faster-report --reads <path_to_reads_folder> [options]

     Required Parameters:
       --reads         Path to the folder containing FASTQ/BAM files.

     Optional Parameters:
       --type          Sequencing platform used ('illumina', 'ont', or 'pacbio').
                       If not set, it is auto-detected from the FASTQ headers.
       --subsample     Fraction of reads to subsample for k-mers calculation (0.1 to 1.0, default: 1.0).
       --outfile       Name of the output HTML report file (default: 'faster-report.html').
       --outdir        Directory where the output report is saved (default: 'output').
       --save_raw      Save raw CSV data used for plotting ('true' or 'false', default: false).

     Metadata Override Options (will override values auto-detected from FASTQ headers):
       --flowcell      Flow cell ID
       --rundate       Run date
       --basecall      Basecaller model
       --user          User name

     Other Options:
       --help          Display this help message.
    ================================================================================
    """
    exit 0
}

if (!params.reads) {
    error "Please provide a path to a folder with fastq/bam files: --reads /path/to/fastq. Use --help for full usage options."
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
    //container 'aangeloo/faster-report'
    //containerOptions "--entrypoint ''"

    input:
    path 'temp/*' // collected there and passed to script

    output:
    path 'header_data.csv'

    script:
    """
    get-header-data.sh temp > header_data.csv
    """
}

process FASTER_REPORT {
    publishDir params.outdir, mode: 'copy'
    container 'aangeloo/faster-report'
    containerOptions "--entrypoint ''"

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
        -r '\\.(fastq|fq|fasta|bam)(\\.gz)?\$' \\
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
