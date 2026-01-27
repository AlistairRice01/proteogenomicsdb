/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAT_FASTQ         } from '../../../modules/nf-core/cat/fastq/main.nf'
include { FASTQC_TRIMGALORE } from '../fastqc_trimgalore/fastqc_trimgalore.nf'
include { MULTIQC           } from '../../../modules/nf-core/multiqc/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN CBIOPORTAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTQC_WORKFLOW {

take:

    //inputs from the rnaseq_workflow   
    samplesheet    //channel: [val(meta), [ samplesheet ] ]
    multiqc_config //channel: /path/to/multiqc config

main:

    //creates empty channels for tool versions, multiqc_files, and the multiqc_report
    multiqc_files_ch  = Channel.empty()
    multiqc_report_ch = Channel.empty()
    versions_ch       = Channel.empty()

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

    //pass the channels into the FASTQC_TRIMGALORE subworkflow
    //this takes the sequencing reads and runs fastqc and trims them
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc,
        params.skip_trimming
    )
    versions_ch       = versions_ch.mix(FASTQC_TRIMGALORE.out.versions_ch)
    multiqc_files_ch  = multiqc_files_ch.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    multiqc_files_ch  = multiqc_files_ch.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        multiqc_files_ch.collect(),
        multiqc_config,
        [],
        [],
        [],
        []
    )
    multiqc_report_ch = MULTIQC.out.report.collect()
    versions_ch       = versions_ch.mix(MULTIQC.out.versions).collect()

emit:

    //emits to rnaseq_workflow 
    multiqc_report_ch //channel: /path/to/multiqc_report.html
    versions_ch       //channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
