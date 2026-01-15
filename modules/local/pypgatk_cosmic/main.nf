process PYPGATK_COSMIC {
    label 'process_medium'
    label 'process_single_thread'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f734a98d2a7b9c84d988474c7e71aef2b0abb987b36fb83b17d52bec0d3babe/data' :
        'community.wave.seqera.io/library/pypgatk_pip:456a1305c1d65d3a' }"

    input:
    path cosmic_config
    val username
    val password

    output:
    path "All_COSMIC_Genes.fasta"    , emit: cosmic_genes
    path "CosmicMutantExport.tsv"    , emit: cosmic_mutations
    path "All_CellLines_Genes.tsv"   , emit: cosmSic_celllines_genes
    path "CosmicCLP_MutantExport.tsv", emit: cosmic_celllines_mutations
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''

    """
    pypgatk_cli.py cosmic-downloader \\
        --config_file ${cosmic_config} \\
        --username ${username} \\
        --password ${password} \\
        --output_directory . \\
 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgatk: \$(echo \$(pypgatk --version 2>&1) | sed 's/^pypgatk v//')
    END_VERSIONS
    """

    stub:
    
    """
    touch All_COSMIC_Genes.fasta
    touch CosmicMutantExport.tsv  
    touch All_CellLines_Genes.tsv
    touch CosmicCLP_MutantExport.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgatk: \$(echo \$(pypgatk --version 2>&1) | sed 's/^pypgatk v//')
    END_VERSIONS
    """
}