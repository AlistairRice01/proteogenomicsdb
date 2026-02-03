/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAT_FASTQ         } from '../../../modules/nf-core/cat/fastq/main.nf'
include { FASTQC as FASTQC_RAW } from '../../../modules/nf-core/fastqc/main.nf'
include { UMITOOLS_EXTRACT  } from '../../../modules/nf-core/umitools/extract/main.nf'
include { TRIMGALORE        } from '../../../modules/nf-core/trimgalore/main.nf'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastq/main.nf'
include { SORTMERNA         } from '../../../modules/nf-core/sortmerna/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN CBIOPORTAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PRE_PROCESSING_WORKFLOW {

take:

    //inputs from the rnaseq_workflow   
    samplesheet    //channel: [val(meta), [ samplesheet ] ]

main:

    //creates empty channels for tool versions, multiqc_files, and the multiqc_report
    multiqc_files_ch  = Channel.empty()
    multiqc_report_ch = Channel.empty()
    fastqc_html       = Channel.empty()
    fastqc_zip        = Channel.empty()
    versions_ch       = Channel.empty()

    trim_reads = Channel.empty()
    trim_html  = Channel.empty()
    trim_zip   = Channel.empty()
    trim_log   = Channel.empty()
    
    //branches the samplesheet into [single] and [multiple] based on how many reads are present 
    samplesheet.branch {
                    meta, fastqs ->
                        single  : fastqs.size() == 1
                            return [ meta, fastqs.flatten() ]
                        multiple: fastqs.size() > 1
                            return [ meta, fastqs.flatten() ]
                }
                .set { ch_fastq }

    //CAT_FASTQ Concatenates FastQ files from same sample if required
    CAT_FASTQ (
        ch_fastq.multiple
    )
    versions_ch = versions_ch.mix(CAT_FASTQ.out.versions_cat).collect()
    ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_fastq.single)

if (!params.skip_fastqc) {

    //pass the channels into the FASTQC_TRIMGALORE subworkflow
    //this takes the sequencing reads and runs fastqc and trims them
    FASTQC_RAW (
        ch_cat_fastq,
    )
    versions_ch = versions_ch.mix(FASTQC_RAW.out.versions_fastqc.first())
    fastqc_html = FASTQC_RAW.out.html
    fastqc_zip  = FASTQC_RAW.out.zip
}

umi_reads = ch_cat_fastq
if (!params.skip_umi) {

    UMITOOLS_EXTRACT (
        ch_cat_fastq
    )
    versions_ch = versions_ch.mix(UMITOOLS_EXTRACT.out.versions.first())
    umi_reads = UMITOOLS_EXTRACT.out.reads
    umi_log   = UMITOOLS_EXTRACT.out.log

}

trim_reads = umi_reads
if (!params.skip_trimming) {

    TRIMGALORE ( 
            umi_reads 
        )
        trim_reads = TRIMGALORE.out.reads
        trim_html  = TRIMGALORE.out.html
        trim_zip   = TRIMGALORE.out.zip
        trim_log   = TRIMGALORE.out.log
        versions_ch   = versions_ch.mix(TRIMGALORE.out.versions_trimgalore.first())


    versions_ch       = versions_ch.mix(FASTQC_TRIMGALORE.out.versions_ch)
    multiqc_files_ch  = multiqc_files_ch.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    multiqc_files_ch  = multiqc_files_ch.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

}

if (!params.skip_fastqc) {

    FASTQC_TRIM (
        trim_reads
    )
    versions_ch = versions_ch.mix(FASTQC_TRIM.out.versions_fastqc.first())
    fastqc_html = FASTQC_TRIM.out.html
    fastqc_zip  = FASTQC_TRIM.out.zip
}

sortmerna_reads = trim_reads
if (!params.skip_sortmerna) {

    SORTMERNA (
        trim_reads,
        [],
        []
    )
    sortmerna_reads = SORTMERNA.out.reads.collect()

}


emit:

    sortmerna_reads   //channel: [ val(meta), [reads] ]
    versions_ch       //channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Read QC, UMI extraction and trimming
//