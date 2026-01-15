
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PYPGATK_COSMIC   } from '../../../modules/local/pypgatk_cosmic/main.nf'
include { PYPGATK_COSMICDB } from '../../../modules/local/pypgatk_cosmicdb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN COSMICDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COSMICDB {

// inputs from the main workflow
take:
    cosmic_config       //channel: contains the cosmic_config file used in both downloading and generating a database from COSMIC data
    username_ch         //channel: contains the users COSMIC username 
    password_ch         //channel: contains the users COSMIC password

main:

    versions_ch = Channel.empty()

    //PYPGATK_COSMIC downlownloads data from the COSMIC database using the cosmic_config, username, and password
    PYPGATK_COSMIC (
        cosmic_config,
        username_ch,
        password_ch
    )
    versions_ch = versions_ch.mix(PYPGATK_COSMIC.out.versions).collect()

    //creates an empty channel that will then be populated by the COSMIC gene data 
    cosmic_genes = Channel.empty()
    cosmic_genes = PYPGATK_COSMIC.out.cosmic_genes.collect()

    //creaes an empty channel that will be populated with the COSMIC mutation data
    comsic_mutations = Channel.empty()
    cosmic_mutations = PYPGATK_COSMIC.out.cosmic_mutations.collect()

    //PYPGATK_COSMICDB generates a database using the COSMIC data downloaded
    PYPGATK_COSMICDB ( 
        cosmic_genes.map { [ [:], it ] },
        cosmic_mutations.map { [ [:], it ] },
        cosmic_config
    )   
    versions_ch = versions_ch.mix(PYPGATK_COSMICDB.out.versions).collect()

    //creates an empty channel that will then be populated with the COSMIC database
    cosmic_database = Channel.empty()
    cosmic_database = PYPGATK_COSMICDB.out.fasta.collect()

// emits to the main workflow
emit:
    cosmic_database     //channel: contains the database generated from this workflow
    versions_ch         //channel: contains versions.yml holding the version information for each of the tools

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/