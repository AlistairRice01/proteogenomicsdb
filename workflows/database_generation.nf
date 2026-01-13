
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC_WORKFLOW } from '../subworkflows/local/fastqc_workflow/fastqc_workflow.nf'

include { RNASEQDB     } from '../subworkflows/local/rnaseq_workflow/rnaseq_workflow.nf'
include { ENSEMBLDB    } from '../subworkflows/local/ensembl_workflow/ensembl_workflow.nf'
include { COSMICDB     } from '../subworkflows/local/cosmic_workflow/cosmic_workflow.nf'
include { GNOMADDB     } from '../subworkflows/local/gnomad_workflow/gnomad_workflow.nf'
include { CBIOPORTALDB } from '../subworkflows/local/cbioportal_workflow/cbioportal_workflow.nf'
include { MERGEDB      } from '../subworkflows/local/merge_workflow/merge_workflow.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DATABASE_GENERATION {

take:
    //proteogenomicsdb
    samplesheet
    reference
    annotation
    transcripts
    custom_config
    dna_config

    //ensembldb
    ensembl_downloader_config
    species_id 
    ensembl_config 
    altorfs_config 
    pseudogenes_config
    ncrna_config 

    //cosmicdb
    cosmic_config
    username
    password

    //gnomad/genecodedb
    genecode_transcripts_url
    genecode_annotations_url
    gnomad_url
    gnomad_config

    //cbioportaldb
    grch37_url
    cbio_config
    cbio_study

    //mergedb
    minimum_aa
    stop_codons
    decoy_config

main:

    //create an empty channel which will later contain the version information for each of the tools
    Channel
        .empty()
        .set { versions }
    
    Channel
        .empty()
        .set { multiqc_report }
    
    FASTQC_WORKFLOW (
    samplesheet,
    versions,
    multiqc_report
    )
    versions = versions.mix(FASTQC_WORKFLOW.out.versions).collect()
    multiqc_report = FASTQC_WORKFLOW.out.multiqc_report.collect() 

    //create an empty channel which will later contain all of the databases mixed together
    Channel
        .empty()
        .set { mixed_databases }


    //conditional execution based on wether the workflow is turned on or off in the config file 
    if (params.proteogenomicsdb) {

        //pass the channels into the PROTEOGENOMICSDB subworkflow - this takes sequencing data to produce a novel protein database
        RNASEQDB (
            samplesheet,
            annotation,
            reference,
            custom_config,
            dna_config,
            transcripts,
            versions
        )
        //extract the version information from the subworkflow
        versions = versions.mix(RNASEQDB.out.versions).collect()
        //extract the databases from PROTEOGENOMICSDB and add them to mixed_databases_ch
        mixed_databases_ch = mixed_databases_ch.mix(RNASEQDB.out.merged_databases).collect()
        //extract the multiqc report from PROTEOGENOMICSDB
        multiqc_report = RNASEQDB.out.multiqc_report.collect()

    } 
    else {
        //bypass the subworkflow
        log.info "proteogenomicsdb subworkflow skipped."
    }

    //conditional execution based on wether the workflow is turned on or off in the config file
    if (params.ensembldb) {
        //pass the channels into the ENSEMBLDB subworkflow - this downloads data FROM ENSEMBL to create a protein database
        ENSEMBLDB (
            ensembl_downloader_config,
            species_id,
            ensembl_config,
            altorfs_config,
            pseudogenes_config,
            ncrna_config,
            versions
        )
        //extract the version information from the subworkflow
        versions = versions.mix(ENSEMBLDB.out.versions).collect()
        //extract the databases from ENSEMBLDB and add them to mixed_databases_ch
        mixed_databases = mixed_databases.mix(ENSEMBLDB.out.mixed_databases).collect()

    }
    else {
        //bypass the subworkflow
        log.info "ensembldb subworkflow skipped."
    }

    //conditional execution based on wether the workflow is turned on or off in the config file
    if (params.cosmicdb) {

        //pass the channels into the COSMICDB subworkflow - this downloads data FROM COSMIC to create a protein database
        COSMICDB (
            cosmic_config,
            username,
            password,
            versions
        )
        //extract the version information from the subworkflow
        versions = versions.mix(COSMICDB.out.versions).collect()   
        //extract the databases from COSMICDB and add them to mixed_databases
        mixed_databases = mixed_databases.mix(COSMICDB.out.cosmic_database).collect()

    }
    else {
        //bypass the subworkflow
        log.info "cosmicdb subworkflow skipped."
    }
    
    //conditional execution based on wether the workflow is turned on or off in the config file
    if (params.gnomaddb) {

        //pass the channels into the GNOMADDB subworkflow - this downloads data FROM GENECODE and GNOMAD to create a protein database
        GNOMADDB (
            genecode_transcripts_url,
            genecode_annotations_url,
            gnomad_url,
            gnomad_config,
            versions
        )
        //extract the version information from the subworkflow
        versions = versions.mix(GNOMADDB.out.versions).collect()
        //extract the databases from GNOMADDB and add them to mixed_databases_ch
        mixed_databases = mixed_databases.mix(GNOMADDB.out.gnomad_database).collect()

    }
    else {
        //bypass the subworkflow
        log.info "gnomaddb subworkflow skipped."
    }

    //conditional execution based on wether the workflow is turned on or off in the config file
    if (params.cbioportaldb) {

        //pass the channels into the CBIOPORTALDB subworkflow - this downloads data FROM CBIOPORTAL to create a protein database
        CBIOPORTALDB (
            grch37_url,
            cbio_config,
            cbio_study,
            versions
        )
        //extract the version information from the subworkflow
        versions = versions.mix(CBIOPORTALDB.out.versions).collect()
        //extract the databases from CBIOPORTALDB and add them to mixed_databases_ch
        mixed_databases_ch = mixed_databases.mix(CBIOPORTALDB.out.cbioportal_database).collect()

    }
    else {
        //bypass the subworkflow
        log.info "cbioportaldb subworkflow skipped."
    }

    //pass the channels into the MERGEDB subworkflow - this merges all of the databases together
    MERGEDB (
        mixed_databases,
        ensembl_config,
        minimum_aa,
        stop_codons,
        decoy_config,
        versions
    )
    //extract the version information from the subworkflow
    versions = versions.mix(MERGEDB.out.versions).collect()

    //extract the databases from MERGEDB and add overwrite mixed_databases
    mixed_databases = MERGEDB.out.databases.collect()
    
    //create a channel containing the decoy database
    Channel
        .empty()
        .set { decoy_database }
    decoy_database = MERGEDB.out.decoy.collect()


emit:

    //emit the final files from the workflow
    mixed_databases      //channel: contains the final peptide database
    decoy_database       //channel: contains the decoy database
    versions             //channel: contains the version information for each of the tools used in the pipeline
    multiqc_report       //channel: contains the multiqc report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

