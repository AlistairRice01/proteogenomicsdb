/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COSMIC_DOWNLOAD           } from '/exports/eddie/scratch/s2215490/test_cosmic/cosmic.nf'
include { PYPGATK_COSMICDB          } from '/exports/eddie/scratch/s2215490/proteogenomicsdb/modules/local/pypgatk/cosmic_to_proteindb/main.nf'
include { PYPGATK_COSMICDB as PYPGATK_COSMICDB_CELLINES } from '/exports/eddie/scratch/s2215490/proteogenomicsdb/modules/local/pypgatk/cosmic_to_proteindb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN COSMICDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    params.cosmic_config = "/exports/eddie/scratch/s2215490/proteogenomicsdb/conf/pypgatk_config_files/cosmic_config.yaml"
    params.username      = "s2215490@ed.ac.uk" 
    params.password      = "BeepBoop8294!"
    params.output_genes  = "Cosmic_Genes_Tsv_v103_GRCh38.tar"
    params.output_mutations  = "Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38.tar"
    params.genes         = "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v103/Cosmic_Genes_Tsv_v103_GRCh38.tar&bucket=downloads"
    params.mutations     = "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cosmic/v103/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38.tar&bucket=downloads"
    params.cosmic_url_celllines_genes = "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cell_lines/v103/CellLinesProject_CompleteGeneExpression_Tsv_v103_GRCh38.tar&bucket=downloads"
    params.cosmic_url_celllines_mutations = "https://cancer.sanger.ac.uk/api/mono/products/v1/downloads/scripted?path=grch38/cell_lines/v103/CellLinesProject_MutationTracking_Tsv_v103_GRCh38.tar&bucket=downloads"
    params.skip_celllines = false
    params.skip_cosmic = false
    params.cosmic_cancer_type = 'all'

workflow {

take:
// inputs from the main workflow

    cosmic_config = Channel.from(params.cosmic_config)
    username_ch = Channel.from(params.username) 
    password_ch = Channel.from(params.password)
    output_file_mutations = Channel.from(params.output_mutations) 
    output_file_genes = Channel.from(params.output_genes)
    cosmic_url_genes = Channel.from(params.genes)
    cosmic_url_mutations = Channel.from(params.mutations)
    cosmic_url_celllines_genes = Channel.from(params.cosmic_url_celllines_genes)
    cosmic_url_celllines_mutations = Channel.from(params.cosmic_url_celllines_mutations)


main:

    versions_ch     = Channel.empty()
    cosmic_database = Channel.empty()

    cosmic_genes        = Channel.empty()
    cosmic_mutations    = Channel.empty()
    celllines_genes     = Channel.empty()
    celllines_mutations = Channel.empty()

    cosmic_genes_untar        = Channel.empty()
    cosmic_mutations_untar    = Channel.empty()
    celllines_genes_untar     = Channel.empty()
    celllines_mutations_untar = Channel.empty()

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

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/