
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WGET   as WGET_TRANSCRIPTS    } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP as GUNZIP_TRANSCRIPTS  } from '../../../modules/nf-core/gunzip/main.nf'
include { WGET as WGET_ANNOTATION       } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP as GUNZIP_ANNOTATION   } from '../../../modules/nf-core/gunzip/main.nf'
include { GSUTIL                        } from '../../../modules/local/gsutil/main.nf'
include { TABIX_BGZIP                   } from '../../../modules/nf-core/tabix/bgzip/main.nf'
include { PYPGATK                       } from '../../../modules/local/pypgatk_vcfdb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN GNOMADDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GNOMADDB {

take:

    //inputs from the main workflow
    genecode_transcripts_url        //channel: contains the url used to download transcripts from genecode
    genecode_annotation_url         //channel: contains the url used to download annotations from genecode
    gnomad_url                      //channel: contains the url used to download a vcf file from gnomad
    ensembl_config                  //channel: contains the ensembl_config used in database generation
    versions_ch                     //channel: contains versions.yml holding the version information for each of the tools

main:

    //WGET_TRANSCRIPTS downloads the transcript files from genecode 
    WGET_TRANSCRIPTS (
        genecode_transcripts_url.map { [ [id: 'transcripts.fa' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_TRANSCRIPTS.out.versions).collect()

    //creates an empty channel which will then be populated with the transcripts downloaded from genecode
    Channel 
        .empty()
        .set { genecode_transcripts_ch }
    genecode_transcripts_ch = WGET_TRANSCRIPTS.out.outfile.collect()

    //GUNZIP_TRANSCRIPTS unzipps the transcript files downloaded from genecode
    GUNZIP_TRANSCRIPTS (
        genecode_transcripts_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_TRANSCRIPTS.out.versions).collect()

    //WGET_ANNOTATION downloads the annotation files from genecode 
    WGET_ANNOTATION (
        genecode_annotation_url.map { [ [id: 'annotations.gtf' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_ANNOTATION.out.versions).collect()

    //creates an empty channel which will then be populated with the annotation downloaded from genecode
    Channel 
        .empty()
        .set { genecode_annotation_ch }
    genecode_annotation_ch = WGET_ANNOTATION.out.outfile.collect()

    //GUNZIP_ANNOTATION unzipps the annotation files downloaded from genecode
    GUNZIP_ANNOTATION (
        genecode_annotation_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_ANNOTATION.out.versions).collect()
    
    //GSUTIL downloads the vcf file from gnomad 
    GSUTIL (
        gnomad_url.map { [ [id: 'gnomad.vcf' ], it ] }
    )
    versions_ch = versions_ch.mix(GSUTIL.out.versions).collect()

    //creates ane mpty channel that will then be populated with the vcf file downloaded from gnomad
    Channel
        .empty()
        .set { vcf_compressed }
    vcf_compressed = GSUTIL.out.vcf.collect()

    //TABIX_BGZIP unzips the vcf file downloaded from gnomad
    TABIX_BGZIP (
        vcf_compressed
    )
    versions_ch = versions_ch.mix(TABIX_BGZIP.out.versions).collect()

    //creates a channel that will then be populated with the unzipped vcf file  
    Channel
        .empty()
        .set { gnomad_vcf_extracted }
    gnomad_vcf_extracted = TABIX_BGZIP.out.output.collect()

    //creates a channel that will then be populated with the unzipped annotation file
    Channel
        .empty()
        .set { genecode_annotation_gunzipped }
    genecode_annotation_gunzipped = GUNZIP_ANNOTATION.out.gunzip.collect()

    //creates a channel that will then be populated with the unzipped transcript file
    Channel
        .empty()
        .set { genecode_transcripts_gunzipped }
    genecode_transcripts_gunzipped = GUNZIP_TRANSCRIPTS.out.gunzip.collect()

    //PYPGATK takes the unzipped vcf, annotation, and transcript along with the ensembl config to generate a peptide database
    PYPGATK (
        gnomad_vcf_extracted,
        genecode_annotation_gunzipped,
        genecode_transcripts_gunzipped,
        ensembl_config  
    )
    versions_ch = versions_ch.mix(PYPGATK.out.versions).collect()

    //creates an empty channel that will be populated with the database generated from PYPGATK
    Channel
        .empty()
        .set { gnomad_database }
    gnomad_database = PYPGATK.out.database.collect()

emit:
    
    // emits to the main workflow
    gnomad_database     //channel: contains the database generated from this workflow
    versions_ch         //channel: contains versions.yml holding the version information for each of the tools

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/