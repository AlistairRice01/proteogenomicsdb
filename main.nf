#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/proteogenomicsdb
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/proteogenomicsdb
    Website: https://nf-co.re/proteogenomicsdb
    Slack  : https://nfcore.slack.com/channels/proteogenomicsdb
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DATABASE_GENERATION     } from './workflows/database_generation.nf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_proteogenomicsdb_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_proteogenomicsdb_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_proteogenomicsdb_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.reference = getGenomeAttribute('fasta')
params.annotation = getGenomeAttribute('gtf')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_PROTEOGENOMICSDB {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // WORKFLOW: Run pipeline
    //
    DATABASE_GENERATION (
        //proteogenomics paramiters
        samplesheet,            //samplesheet containing the rna-seq transcripts
        params.reference,       //
        params.annotation,      //
        params.transcripts,     //
        params.custom_config,   //path to the custom_config
        params.dna_config,      //path to the dna_config
        //ensembl paramaters
        params.ensembl_downloader_config,   //path to the ensembl_downloader_config
        params.species_id,                  //species id to download from ensembl
        params.ensembl_config,              //path to the ensembl_config
        params.altorfs_config,              //path to the altorfs_config
        params.pseudogenes_config,          //path to the pseudogenes_config
        params.ncrna_config,                //path to the ncrna_config
        //cosmic paramiters
        params.cosmic_config,   //path to the cosmic_config
        params.username,        //            
        params.password,        //
        //genecode/gnomad paramaters
        params.genecode_transcripts_url,    //
        params.genecode_annotations_url,    //
        params.gnomad_url,                  //
        params.gnomad_config,               //path to the gnomad_config
        //cbioportal paramaters
        params.grch37_url,      //      
        params.cbio_config,     //path to the cbioportal_config
        params.cbio_study,      //
        //merge paramaters
        params.minimum_aa,      //
        params.stop_codons,     //
        params.decoy_config     //path to the decoy_config

    )
    emit:
    mixed_databases = DATABASE_GENERATION.out.mixed_databases    //channel: contains the final peptide database
    decoy_database  = DATABASE_GENERATION.out.decoy_database     //channel: contains the decoy database
    versions        = DATABASE_GENERATION.out.versions           //channel: contains the version information for each of the tools used in the pipeline
    multiqc_report  = DATABASE_GENERATION.out.multiqc_report     // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_PROTEOGENOMICSDB (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_PROTEOGENOMICSDB.out.multiqc_report,
        NFCORE_PROTEOGENOMICSDB.out.mixed_databases,  
        NFCORE_PROTEOGENOMICSDB.out.decoy_database,     
        NFCORE_PROTEOGENOMICSDB.out.versions 
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
