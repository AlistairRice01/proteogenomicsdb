//
// Read QC, UMI extraction and trimming
//

include { FASTQC           } from '../../modules/nf-core/fastqc/main.nf'
include { TRIMGALORE       } from '../../modules/local/trimgalore/main.nf'

workflow FASTQC_TRIMGALORE {
    
take:
   
    reads         // channel: [ val(meta), [ reads ] ]
    skip_fastqc   // boolean: true/false
    skip_trimming // boolean: true/false
    versions

main:

    Channel 
        .empty()
        .set { fastqc_html }
    
    Channel
        .empty()
        .set { fastqc_zip }
    
    if (!skip_fastqc) {
        FASTQC ( 
            reads.view() 
        )
        versions = versions.mix(FASTQC.out.versions.first())
        fastqc_html = FASTQC.out.html
        fastqc_zip  = FASTQC.out.zip
    }

    Channel 
        .empty()
        .set { trim_reads }
    trim_reads = reads

    Channel
        .empty()
        .set { trim_html }

    Channel
        .empty()
        .set { trim_zip }

    Channel
        .empty()
        .set { trim_log }

    if (!skip_trimming) {
        TRIMGALORE ( 
            reads.view() 
        )
        trim_reads = TRIMGALORE.out.reads
        trim_html  = TRIMGALORE.out.html
        trim_zip   = TRIMGALORE.out.zip
        trim_log   = TRIMGALORE.out.log
        versions   = versions.mix(TRIMGALORE.out.versions_trimgalore.first())
    }

emit:

    reads = trim_reads     // channel: [ val(meta), [ reads ] ]

    fastqc_html            // channel: [ val(meta), [ html ] ]
    fastqc_zip             // channel: [ val(meta), [ zip ] ]

    trim_html              // channel: [ val(meta), [ html ] ]
    trim_zip               // channel: [ val(meta), [ zip ] ]
    trim_log               // channel: [ val(meta), [ txt ] ]

    versions // channel: [ versions.yml ]
}