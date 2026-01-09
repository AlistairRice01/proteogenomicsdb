
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WGET as WGET_GRCH37         } from '../../../modules/nf-core/wget/main.nf'
include { GUNZIP                      } from '../../../modules/nf-core/gunzip/main.nf'
include { PYPGATK_CBIOPORTAL_DOWNLOAD } from '../../../modules/local/pypgatk_download_cbioportal/main.nf'
include { PYPGATK_CBIOPORTAL          } from '../../../modules/local/pypgatk_cbioportaldb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN CBIOPORTAL WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CBIOPORTALDB {

take:

    //inputs from the main workflow
    grch37_link         //channel: html link to the file that you want to download in the grch37 database 
    cbioportal_config   //channel: config file for both downloading and creating a database using cBioPortal data 
    sample_id           //channel: string telling the program which samples to download 
    versions_ch         //channel: contains versions.yml holding the version information for each of the tools

main:

    //WGET_GRCH37 takes the link and downloads it off the internet
     WGET_GRCH37 (
        grch37_link.map { [ [id: 'grch37' ], it ] }
    )
    versions_ch = versions_ch.mix(WGET_GRCH37.out.versions).collect()

    //PYPGATK_CBIOPORTAL_DOWNLOAD downloads the specified samples fromt he cBioPortal database
    PYPGATK_CBIOPORTAL_DOWNLOAD (
        cbioportal_config.map { [ [id: 'cbioportal_config'], it ] },
        sample_id
    )
    versions_ch = versions_ch.mix(PYPGATK_CBIOPORTAL_DOWNLOAD.out.versions).collect()

    //creates an empty channel that will then be populated with the file downloaded off the grch37 database 
    Channel
        .empty()
        .set { grch37_ch }
    grch37_ch = WGET_GRCH37.out.outfile.collect()

    //GUNZIP unzips the file downloaded using WGET_GRCH37 
    GUNZIP(
        grch37_ch
    )
    versions_ch = versions_ch.mix(GUNZIP.out.versions).collect()

    //creates an empty channel which will be populated with the unzipped grch37 file 
    Channel
        .empty()
        .set { grch37_unzipped }
    grch37_unzipped = GUNZIP.out.gunzip.collect()
   
    //creates an empty channel that will then be populated with the cBioPortal mutations file 
    Channel
        .empty()
        .set { cbio_mutations }
    cbio_mutations = PYPGATK_CBIOPORTAL_DOWNLOAD.out.cbio_mutations.collect()

    //creates an empty channel that will then be populated with the cBioPortal samples file 
    Channel
        .empty()
        .set { cbio_samples }
    cbio_samples = PYPGATK_CBIOPORTAL_DOWNLOAD.out.cbio_samples.collect()

    //PYPGATK_CBIOPORTAL using the grch37, mutations, and samples file generates a peptide database 
    PYPGATK_CBIOPORTAL (
        grch37_unzipped,
        cbio_mutations.map { [ [id: 'cbio_mutations' ], it ] },
        cbio_samples.map { [ [id: 'cbio_samples' ], it ] },
        cbioportal_config
    )
    versions_ch = versions_ch.mix(PYPGATK_CBIOPORTAL.out.versions).collect()

    //creates an empty channel that will be populated with the database generated from PYPGATK_CBIOPORTAL
    Channel
        .empty()
        .set { cbioportal_database }
    cbioportal_database = PYPGATK_CBIOPORTAL.out.cbioportal_proteindb.collect()

emit:
    
    //emits to the main workflow 
    cbioportal_database     //channel: contains the database generated from this workflow
    versions_ch             //channel: contains versions.yml holding the version information for each of the tools
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/