//
// Read QC, UMI extraction and trimming
//

include { FASTQC           } from '../../modules/nf-core/fastqc/main.nf'
include { TRIMGALORE       } from '../../modules/nf-core/trimgalore/main.nf'

workflow FASTQC_TRIMGALORE {
    
take:
   
    reads         // channel: [ val(meta), [ reads ] ]
    skip_fastqc   // boolean: true/false
    skip_trimming // boolean: true/false

main:

    fastqc_html = Channel.empty()
    fastqc_zip  = Channel.empty()
    versions_ch = Channel.empty()
    
    if (!skip_fastqc) {
    
        FASTQC ( 
            reads 
        )
        versions_ch = versions_ch.mix(FASTQC.out.versions_fastqc.first())
        fastqc_html = FASTQC.out.html
        fastqc_zip  = FASTQC.out.zip
    }

    else {
        //bypass the subworkflow
        log.info "fastqc skipped."
    }

    trim_reads = Channel.empty()
    trim_html  = Channel.empty()
    trim_zip   = Channel.empty()
    trim_log   = Channel.empty()

    if (!skip_trimming) {
        TRIMGALORE ( 
            reads 
        )
        trim_reads = TRIMGALORE.out.reads
        trim_html  = TRIMGALORE.out.html
        trim_zip   = TRIMGALORE.out.zip
        trim_log   = TRIMGALORE.out.log
        versions_ch   = versions_ch.mix(TRIMGALORE.out.versions_trimgalore.first())
    }

    else {
        //bypass the subworkflow
        log.info "trimgalore skipped."
    }

emit:

    reads = trim_reads     // channel: [ val(meta), [ reads ] ]

    fastqc_html            // channel: [ val(meta), [ html ] ]
    fastqc_zip             // channel: [ val(meta), [ zip ] ]

    trim_html              // channel: [ val(meta), [ html ] ]
    trim_zip               // channel: [ val(meta), [ zip ] ]
    trim_log               // channel: [ val(meta), [ txt ] ]

    versions_ch // channel: [ versions.yml ]
}