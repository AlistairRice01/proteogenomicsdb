include { CAT_FASTQ         } from '../../../modules/nf-core/cat/fastq/main.nf'
include { FASTQC_TRIMGALORE } from '../../nf-core/fastqc_trimgalore.nf'
include { MULTIQC           } from '../../../modules/nf-core/multiqc/main.nf'

workflow FASTQC_WORKFLOW {

take:

    samplesheet
    multiqc_config

main:

    multiqc_files_ch  = Channel.empty()
    multiqc_report_ch = Channel.empty()
    versions_ch       = Channel.empty()

    samplesheet.branch {
                    meta, fastqs ->
                        single  : fastqs.size() == 1
                            return [ meta, fastqs.flatten() ]
                        multiple: fastqs.size() > 1
                            return [ meta, fastqs.flatten() ]
                }
                .set { ch_fastq }

    // MODULE:
    // Concatenate FastQ files from same sample if required
    CAT_FASTQ (
        ch_fastq.multiple
    )
    versions_ch = versions_ch.mix(CAT_FASTQ.out.versions_cat).collect()
    ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_fastq.single)



    // MODULE: Run FastQC, trimgalore!
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc,
        params.skip_trimming,
        versions_ch
    )
    versions_ch       = versions_ch.mix(FASTQC_TRIMGALORE.out.versions_ch)
    multiqc_files_ch  = multiqc_files_ch.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    multiqc_files_ch  = multiqc_files_ch.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    multiqc_config_ch = Channel.fromPath(multiqc_config)

    MULTIQC (
        multiqc_files_ch.collect(),
        multiqc_config_ch,
        [],
        [],
        [],
        []
    )
    multiqc_report_ch = MULTIQC.out.report.collect()
    versions_ch       = versions_ch.mix(MULTIQC.out.versions).collect()

emit:

    multiqc_report_ch
    versions_ch 


}
