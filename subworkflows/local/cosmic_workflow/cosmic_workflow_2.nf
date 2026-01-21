/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COSMIC_DOWNLOAD           } from "$projectDir/modules/local/pypgatk/cosmic_downloader/new_cosmic.nf"
include { PYPGATK_COSMICDB          } from "$projectDir/modules/local/pypgatk/cosmic_to_proteindb/main.nf"
include { PYPGATK_COSMICDB as PYPGATK_COSMICDB_CELLINES } from "$projectDir/modules/local/pypgatk/cosmic_to_proteindb/main.nf"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN COSMICDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COSMICDB {

take:

    cosmic_config
    username_ch
    password_ch
    cosmic_url_genes
    cosmic_url_mutations
    cosmic_url_celllines_genes
    cosmic_url_celllines_mutations

main:

// inputs from the main workflow

    cosmic_config = Channel.from(params.cosmic_config)
    username_ch = Channel.from(params.username) 
    password_ch = Channel.from(params.password)
    output_file_mutations = Channel.from(params.output_mutations) 
    output_file_genes = Channel.from(params.output_genes)
    cosmic_url_genes = Channel.from(params.cosmic_genes_url)
    cosmic_url_mutations = Channel.from(params.cosmic_mutations_url)
    cosmic_url_celllines_genes = Channel.from(params.cosmic_celllines_genes_url)
    cosmic_url_celllines_mutations = Channel.from(params.cosmic_celllines_mutations_url)

    versions_ch = Channel.empty()
    cosmic_database = Channel.empty()

    cosmic_genes = Channel.empty()
    cosmic_mutations = Channel.empty()
    celllines_genes = Channel.empty()
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
        cosmic_genes.map { [ [id: 'genes'], it ] },
        cosmic_mutations.map { [ [id: 'mutations'], it ] },
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
        celllines_genes.map { [ [id: 'cellline_genes'], it ] },
        celllines_mutations.map { [ [id: 'cellline_mutations'], it ] },
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

    cosmic_database     //channel: contains the database generated from this workflow
    versions_ch         //channel: contains versions.yml holding the version information for each of the tools

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/