include { CAT_FASTQ         } from '../../../modules/nf-core/cat/fastq/main.nf'
include { FASTQC_TRIMGALORE } from '../../nf-core/fastqc_trimgalore.nf'
include { MULTIQC           } from '../../../modules/nf-core/multiqc/main.nf'

workflow FASTQC_WORKFLOW {

take:

    samplesheet
    versions
    multiqc_report

main:

    Channel
        .empty()
        .set { ch_multiqc_files }

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
    versions = versions.mix(CAT_FASTQ.out.versions_cat).collect()
    ch_cat_fastq = CAT_FASTQ.out.reads.mix(ch_fastq.single)

    Channel 
        .from(params.skip_fastqc)
        .set { skip_fastqc }

    Channel 
        .from(params.skip_trimming)
        .set { skip_trimming }

        // MODULE: Run FastQC, trimgalore!
    FASTQC_TRIMGALORE (
        ch_cat_fastq.view(),
        skip_fastqc,
        skip_trimming,
        versions
    )
    versions = versions.mix(FASTQC_TRIMGALORE.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files  = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        [],
        [],
        [],
        [],
        []
    )
    multiqc_report = MULTIQC.out.report.collect()
    versions = versions.mix(MULTIQC.out.versions).collect()

emit:

    multiqc_report
    versions 


}
