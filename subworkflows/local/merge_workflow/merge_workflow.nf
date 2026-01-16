
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAT_CAT       } from '../../../modules/nf-core/cat/cat/main.nf'
include { PYPGATK_CLEAN } from '../../../modules/local/pypgatk_clean/main.nf'
include { PYPGATK_DECOY } from '../../../modules/local/pypgatk_decoy/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN MERGEDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MERGEDB {

take:

    // inputs from the main workflow
    mixed_databases     //channel: contains all of the databases generated from the other workflows mixed into a single channel
    clean_config        //channel: contains contains a config file defining how the database should be cleaned
    decoy_config        //channel: contains constains a config file defining how the decoy database should be generated 

main:

    versions_ch = Channel.empty()
    decoy       = Channel.empty()

    //creates an empty channel that will be populated with the databases generated in other workflows and flattened into a list - this is to remove the metadata from each database
    flat_databases = Channel.empty()
    flat_databases = mixed_databases.flatten().filter{it instanceof java.nio.file.Path && it.isFile()}

    //creates an empty channel that will take the flattened databases and collect them 
    collected_databases = Channel.empty()
    collected_databases = flat_databases.collect()

    //CAT_CAT concatenates all of the databases into a single database
    CAT_CAT (
        collected_databases.map { [ [id: 'final_database'], it ] }
    )
    versions_ch = versions_ch.mix(CAT_CAT.out.versions_cat).collect()

    //creates an empty channel that will then be populated with the concatenated database
    databases = Channel.empty()
    databases = CAT_CAT.out.file_out.collect()
        .map { meta, it ->
            return [it] }

    //PYPGATK_CLEAN takes the concatenated database and cleans it
    PYPGATK_CLEAN (
        databases.map { [ [id: 'cleaned_database'], it ] },
        clean_config
    )
    databases = PYPGATK_CLEAN.out.clean_database.collect()
    versions_ch = versions_ch.mix(PYPGATK_CLEAN.out.versions).collect()

if (!params.skip_decoy) {

    //PYPGATK_DECOY generates a decoy database from the cleaned database using the decoy_config
    PYPGATK_DECOY (
        databases.map { [ [id: 'decoy_database'], it ] },
        decoy_config
    )
    versions_ch = versions_ch.mix(PYPGATK_DECOY.out.versions).collect()
    decoy       = PYPGATK_DECOY.out.decoy_database.collect()


}

else {
        //bypass the subworkflow
        log.info "decoy generation skipped."
    }
    //creates an empty channel that will then be populated with the decoy database 

emit:

    // emits to the main workflow
    databases    //channel: contains the final database 
    decoy           //channel: contains the decoy database 
    versions_ch     //channel: contains versions.yml holding the version information for each of the tools 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/