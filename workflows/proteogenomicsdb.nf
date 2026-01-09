
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROTEOGENOMICSDB } from '../subworkflows/local/proteogenomics_workflow/proteogenomics_workflow.nf'
include { ENSEMBLDB        } from '../subworkflows/local/ensembl_workflow/ensembl_workflow.nf'
include { COSMICDB         } from '../subworkflows/local/cosmic_workflow/cosmic_workflow.nf'
include { GNOMADDB         } from '../subworkflows/local/gnomad_workflow/gnomad_workflow.nf'
include { CBIOPORTALDB     } from '../subworkflows/local/cbioportal_workflow/cbioportal_workflow.nf'
include { MERGEDB          } from '../subworkflows/local/merge_workflow/merge_workflow.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROTEOGENOMICS {

main:
    
    //create an empty channel which will later contain all of the databases mixed together
    Channel
        .empty()
        .set { mixed_databases_ch }

    //create an empty channel which will later contain the version information for each of the tools
    Channel
        .empty()
        .set { versions_ch }

    //create a channel containing the PROTEOGENOMICSDB read data
    Channel
        .fromPath(params.reads)
        .set { reads_ch }

    //create a channel containing the PROTEOGENOMICSDB annotation data
    Channel
        .fromPath(params.gtf)
        .set { gtf_ch }
    
    //create a channel containing the PROTEOGENOMICSDB reference data
    Channel
        .fromPath(params.reference)
        .set { reference_ch }
    
    //create a channel containing the PROTEOGENOMICSDB transcript data
    channel 
        .fromPath(params.cdna)
        .set { cdna_ch }

    //create a channel containing a custom config for variant database generation
    Channel
        .fromPath(params.custom_config)
        .set { custom_config_ch }

    //create a channel containing a config file for generating a database using transcript sequences 
    Channel
        .fromPath(params.dna_config)
        .set { dna_config_ch }

    //conditional execution based on wether the workflow is turned on or off in the config file 
    if (params.proteogenomicsdb) {

        //pass the channels into the PROTEOGENOMICSDB subworkflow - this takes sequencing data to produce a novel protein database
        PROTEOGENOMICSDB (
            reads_ch,
            gtf_ch,
            reference_ch,
            custom_config_ch,
            dna_config_ch,
            cdna_ch,
            versions_ch
        )
        //extract the version information from the subworkflow
        versions_ch = versions_ch.mix(PROTEOGENOMICSDB.out.versions_ch).collect()
        //extract the databases from PROTEOGENOMICSDB and add them to mixed_databases_ch
        mixed_databases_ch = mixed_databases_ch.mix(PROTEOGENOMICSDB.out.merged_databases_ch).collect()

    } 
    else {
        //bypass the subworkflow
        log.info "proteogenomicsdb subworkflow skipped."
    }


    //create a channel containing the ensembl downloader config
    Channel
        .fromPath(params.ensembl_downloader)
        .set { ensembl_downloader_ch }

    //create a channel containing the species identifier to download off of ensembl 
    Channel
        .from(params.species)
        .set { species_ch }
    
    //create a channel containg the ensembl config which is used in generating a database from the downloaded ensembl database
    Channel
        .fromPath(params.ensembl)
        .set { ensembl_config_ch }

    //create a channel containing the altorfs config which is used in generating a database from the downloaded altorf files
    Channel
        .fromPath(params.altorfs)
        .set { altorfs_config_ch }

    //create a channel containing the altorfs config which is used in generating a database from the downloaded pseudogenes files
    Channel
        .fromPath(params.pseudogenes)
        .set { pseudogenes_config_ch }
    
    //create a channel containing the altorfs config which is used in generating a database from the downloaded ncrna files
    Channel
        .fromPath(params.ncrna)
        .set { ncrna_config_ch }

    //conditional execution based on wether the workflow is turned on or off in the config file 
    if (params.ensembldb) {

        //pass the channels into the ENSEMBLDB subworkflow - this downloads data FROM ENSEMBL to create a protein database
        ENSEMBLDB (
            ensembl_downloader_ch,
            species_ch,
            ensembl_config_ch,
            altorfs_config_ch,
            pseudogenes_config_ch,
            ncrna_config_ch,
            versions_ch
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

    //create a channel containing the cosmic config which both helps in installation of COSMIC data and building a database from the installed data
    Channel
        .fromPath(params.cosmic)
        .set { cosmic_config_ch }
    
    //create a channel containing the users cosmic username 
    Channel
        .from(params.username)
        .set { username_ch }

    //create a channel containing the users cosmic password
    Channel
        .from(params.password)
        .set { password_ch }

    //conditional execution based on wether the workflow is turned on or off in the config file
    if (params.cosmicdb) {

        //pass the channels into the COSMICDB subworkflow - this downloads data FROM COSMIC to create a protein database
        COSMICDB (
            cosmic_config_ch,
            username_ch,
            password_ch,
            versions_ch
        )
        //extract the version information from the subworkflow
        versions_ch = versions_ch.mix(COSMICDB.out.versions_ch).collect()   
        //extract the databases from COSMICDB and add them to mixed_databases_ch
        mixed_databases_ch = mixed_databases_ch.mix(COSMICDB.out.cosmic_database).collect()

    }
    else {
        //bypass the subworkflow
        log.info "cosmicdb subworkflow skipped."
    }

    //create a channel containing the genecode transcripts url   
    Channel
        .from(params.genecode_transcripts)
        .set { genecode_transcripts_url }

    //create a channel containing the genecode annotations url
    Channel
        .from(params.genecode_annotation)
        .set { genecode_annotations_url }

    //create a channel containing the gnomad url
    Channel
        .from(params.gnomad_vcf)
        .set { gnomad_url }

    //create a channel containing the genomad cofig which is involved in building a database from the downloaded data
    Channel
        .fromPath(params.gnomad_config)
        .set { gnomad_config_ch }
    
    //conditional execution based on wether the workflow is turned on or off in the config file
    if (params.gnomaddb) {

        //pass the channels into the GNOMADDB subworkflow - this downloads data FROM GENECODE and GNOMAD to create a protein database
        GNOMADDB (
            genecode_transcripts_url,
            genecode_annotations_url,
            gnomad_url,
            gnomad_config_ch,
            versions_ch
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


    //create a channel containing the grch37 url
    Channel
        .from(params.grch37)
        .set { grch37_ch }
    
    //create a channel containing the cbioportal config file
    Channel
        .fromPath(params.cbio_config)
        .set { cbio_config_ch }
    
    //creates a channel containing the information on which study to download from cBioPortal 
    Channel
        .from(params.studies)
        .set { cbio_study_ch }

    //conditional execution based on wether the workflow is turned on or off in the config file
    if (params.cbioportaldb) {

        //pass the channels into the CBIOPORTALDB subworkflow - this downloads data FROM CBIOPORTAL to create a protein database
        CBIOPORTALDB (
            grch37_ch,
            cbio_config_ch,
            cbio_study_ch,
            versions_ch
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

    //create a channel containing the minimum number of amino acids which should be considered a protein
    Channel
        .from(params.minimum)
        .set { minimum_aa }
    
    //create a channel containing a boolean statment - tells the tool to add a new protein into the database whenever a stop codon is found
    Channel
        .from(params.stop_codons)
        .set { stop_codons }

    //create a channel containing the decoy config used in the generation of a decoy database
    Channel
        .fromPath(params.decoy)
        .set { decoy_config_ch }

    //pass the channels into the MERGEDB subworkflow - this merges all of the databases together
    MERGEDB (
        mixed_databases_ch,
        ensembl_config_ch,
        minimum_aa,
        stop_codons,
        decoy_config_ch,
        versions_ch
    )
    //extract the version information from the subworkflow
    versions_ch = versions_ch.mix(MERGEDB.out.versions_ch).collect()

    //extract the databases from MERGEDB and add overwrite mixed_databases_ch
    mixed_databases_ch = MERGEDB.out.databases_ch.collect()
    
    //create a channel containing the decoy database
    Channel
        .empty()
        .set { decoy_database_ch }
    decoy_database_ch = MERGEDB.out.decoy.collect()


emit:

    //emit the final files from the workflow
    mixed_databases_ch      //channel: contains the final peptide database
    decoy_database_ch       //channel: contains the decoy database
    versions_ch             //channel: contains the version information for each of the tools used in the pipeline

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

