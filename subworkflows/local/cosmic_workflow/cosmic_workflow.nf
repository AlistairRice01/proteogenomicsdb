/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COSMIC_DOWNLOAD           } from '../../../modules/local/cosmic_downloader/main.nf'
include { PYPGATK_COSMICDB          } from '../../../modules/local/pypgatk/cosmic_to_proteindb/main.nf'
include { PYPGATK_COSMICDB as PYPGATK_COSMICDB_CELLINES } from '../../../modules/local/pypgatk/cosmic_to_proteindb/main.nf'
'../../../modules/nf-core/wget/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN COSMICDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COSMICDB {

take:
// inputs from the main workflow

    cosmic_config
    username_ch
    password_ch
    cosmic_url_genes
    cosmic_url_mutations
    cosmic_url_celllines_genes
    cosmic_url_celllines_mutations


main:

    versions_ch     = Channel.empty()
    cosmic_database = Channel.empty()

    cosmic_genes        = Channel.empty()
    cosmic_mutations    = Channel.empty()
    celllines_genes     = Channel.empty()
    celllines_mutations = Channel.empty()

    //PYPGATK_COSMIC downlownloads data from the COSMIC database using the cosmic_config, username, and password
    COSMIC_DOWNLOAD (
        username_ch,
        password_ch,
        cosmic_url_genes,
        cosmic_url_mutations,
        cosmic_url_celllines_genes,
        cosmic_url_celllines_mutations
    )
    //creates an empty channel that will then be populated by the COSMIC gene data 
    cosmic_genes        = COSMIC_DOWNLOAD.out.cosmic_genes.collect()
    cosmic_mutations    = COSMIC_DOWNLOAD.out.cosmic_mutants.collect()
    celllines_genes     = COSMIC_DOWNLOAD.out.cosmic_celllines_genes.collect()
    celllines_mutations = COSMIC_DOWNLOAD.out.cosmic_celllines_mutants.collect()

    if (!params.skip_cosmic) {

    //PYPGATK_COSMICDB generates a database using the COSMIC data downloaded
    PYPGATK_COSMICDB ( 
        cosmic_genes.map { [ [:], it ] },
        cosmic_mutations.map { [ [:], it ] },
        cosmic_config
    )   
    versions_ch = versions_ch.mix(PYPGATK_COSMICDB.out.versions).collect()
    cosmic_database = PYPGATK_COSMICDB.out.database.collect()

    }

    else {
        //bypass the subworkflow
        log.info "cosmic database skipped."
    }

    if (!params.skip_celllines) {

    PYPGATK_COSMICDB_CELLINES ( 
        celllines_genes.map { [ [:], it ] },
        celllines_mutations.map { [ [:], it ] },
        cosmic_config
    )   
    versions_ch = versions_ch.mix(PYPGATK_COSMICDB.out.versions).collect()
    cosmic_database = cosmic_database.mix(PYPGATK_COSMICDB.out.database).collect()

    }

    else {
        //bypass the subworkflow
        log.info "cosmic celllines database skipped."
    }

emit:

    cosmic_database
    versions_ch

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/