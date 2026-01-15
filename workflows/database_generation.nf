
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQDB        } from '../subworkflows/local/rnaseq_workflow/rnaseq_workflow.nf'
include { ENSEMBLDB       } from '../subworkflows/local/ensembl_workflow/ensembl_workflow.nf'
include { COSMICDB        } from '../subworkflows/local/cosmic_workflow/cosmic_workflow.nf'
include { GNOMADDB        } from '../subworkflows/local/gnomad_workflow/gnomad_workflow.nf'
include { CBIOPORTALDB    } from '../subworkflows/local/cbioportal_workflow/cbioportal_workflow.nf'
include { MERGEDB         } from '../subworkflows/local/merge_workflow/merge_workflow.nf'

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
    clean_config
    decoy_config

main:

    //create an empty channel which will later contain the version information for each of the tools
    versions_ch        = Channel.empty()
    mixed_databases_ch = Channel.empty()
    multiqc_report_ch  = channel.empty()

    //conditional execution based on wether the workflow is turned on or off in the config file 
    if (!params.skip_rnaseqdb) {

        reference_ch     = Channel.fromPath(reference)
        annotation_ch    = Channel.fromPath(annotation)
        transcripts_ch   = Channel.fromPath(transcripts)
        custom_config_ch = Channel.fromPath(custom_config)
        dna_config_ch    = Channel.fromPath(dna_config)

        //pass the channels into the PROTEOGENOMICSDB subworkflow - this takes sequencing data to produce a novel protein database
        RNASEQDB (
            samplesheet,
            annotation_ch,
            reference_ch,
            custom_config_ch,
            dna_config_ch,
            transcripts_ch
        )
        //extract the version information from the subworkflow
        versions_ch        = versions_ch.mix(RNASEQDB.out.versions_ch).collect()
        mixed_databases_ch = mixed_databases_ch.mix(RNASEQDB.out.merged_databases_ch).collect()
        multiqc_report_ch  = RNASEQDB.out.multiqc_report_ch

    } 
    else {
        //bypass the subworkflow
        log.info "proteogenomicsdb subworkflow skipped."
    }

    //conditional execution based on wether the workflow is turned on or off in the config file
    if (!params.skip_ensembldb) {

        ensembl_downloader_config_ch = Channel.fromPath(ensembl_downloader_config)        
        ensembl_config_ch            = Channel.fromPath(ensembl_config)    
        altorfs_config_ch            = Channel.fromPath(altorfs_config)
        pseudogenes_config_ch        = Channel.fromPath(pseudogenes_config)
        ncrna_config_ch              = Channel.fromPath(ncrna_config)
        species_id_ch                = Channel.from(species_id)

        //pass the channels into the ENSEMBLDB subworkflow - this downloads data FROM ENSEMBL to create a protein database
        ENSEMBLDB (
            ensembl_downloader_config_ch,
            species_id_ch,
            ensembl_config_ch,
            altorfs_config_ch,
            pseudogenes_config_ch,
            ncrna_config_ch
        )
        //extract the version information from the subworkflow
        versions_ch = versions_ch.mix(ENSEMBLDB.out.versions_ch).collect()
        //extract the databases from ENSEMBLDB and add them to mixed_databases_ch
        mixed_databases_ch = mixed_databases_ch.mix(ENSEMBLDB.out.mixed_databases).collect()

    }
    else {
        //bypass the subworkflow
        log.info "ensembldb subworkflow skipped."
    }

    //conditional execution based on wether the workflow is turned on or off in the config file
    if (!params.skip_cosmicdb) {

        cosmic_config_ch = Channel.fromPath(cosmic_config)
        username_ch      = Channel.from(username)
        password_ch      = Channel.from(password)

        //pass the channels into the COSMICDB subworkflow - this downloads data FROM COSMIC to create a protein database
        COSMICDB (
            cosmic_config_ch,
            username_ch,
            password_ch
        )
        //extract the version information from the subworkflow
        versions_ch = versions_ch.mix(COSMICDB.out.versions_ch).collect()   
        //extract the databases from COSMICDB and add them to mixed_databases
        mixed_databases_ch = mixed_databases_ch.mix(COSMICDB.out.cosmic_database).collect()

    }
    else {
        //bypass the subworkflow
        log.info "cosmicdb subworkflow skipped."
    }
    
    //conditional execution based on wether the workflow is turned on or off in the config file
    if (!params.skip_gnomaddb) {

        genecode_transcripts_url_ch = Channel.from(genecode_transcripts_url)
        genecode_annotations_url_ch = Channel.from(genecode_annotations_url)        
        gnomad_url_ch               = Channel.from(gnomad_url)
        gnomad_config_ch            = Channel.fromPath(gnomad_config)
        
        //pass the channels into the GNOMADDB subworkflow - this downloads data FROM GENECODE and GNOMAD to create a protein database
        GNOMADDB (
            genecode_transcripts_url_ch,
            genecode_annotations_url_ch,
            gnomad_url_ch,
            gnomad_config_ch
        )
        //extract the version information from the subworkflow
        versions_ch = versions_ch.mix(GNOMADDB.out.versions_ch).collect()
        //extract the databases from GNOMADDB and add them to mixed_databases_ch
        mixed_databases_ch = mixed_databases_ch.mix(GNOMADDB.out.gnomad_database).collect()

    }
    else {
        //bypass the subworkflow
        log.info "gnomaddb subworkflow skipped."
    }

    //conditional execution based on wether the workflow is turned on or off in the config file
    if (!params.skip_cbioportaldb) {

        grch37_url_ch  = Channel.from(grch37_url)
        cbio_study_ch  = Channel.from(cbio_study)
        cbio_config_ch = Channel.fromPath(cbio_config)
        
        //pass the channels into the CBIOPORTALDB subworkflow - this downloads data FROM CBIOPORTAL to create a protein database
        CBIOPORTALDB (
            grch37_url_ch,
            cbio_config_ch,
            cbio_study_ch
        )
        //extract the version information from the subworkflow
        versions_ch = versions_ch.mix(CBIOPORTALDB.out.versions_ch).collect()
        //extract the databases from CBIOPORTALDB and add them to mixed_databases_ch
        mixed_databases_ch = mixed_databases_ch.mix(CBIOPORTALDB.out.cbioportal_database).collect()

    }
    else {
        //bypass the subworkflow
        log.info "cbioportaldb subworkflow skipped."
    }

        minimum_aa_ch     = Channel.from(minimum_aa)
        stop_codons_ch    = Channel.from(stop_codons)
        clean_config_ch = Channel.from(clean_config)
        decoy_config_ch   = Channel.fromPath(decoy_config)

    //pass the channels into the MERGEDB subworkflow - this merges all of the databases together
    MERGEDB (
        mixed_databases_ch,
        clean_config_ch,
        minimum_aa_ch,
        stop_codons_ch,
        decoy_config_ch
    )
    //extract the version information from the subworkflow
    versions_ch = versions_ch.mix(MERGEDB.out.versions_ch).collect()

    //extract the databases from MERGEDB and add overwrite mixed_databases
    mixed_databases_ch = MERGEDB.out.databases.collect()
    
    //create a channel containing the decoy database
    decoy_database_ch = Channel.empty()
    decoy_database_ch = MERGEDB.out.decoy.collect()


emit:

    //emit the final files from the workflow
    mixed_databases_ch   //channel: contains the final peptide database
    decoy_database_ch    //channel: contains the decoy database
    versions_ch          //channel: contains the version information for each of the tools used in the pipeline
    multiqc_report_ch    //channel: contains the multiqc report

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

