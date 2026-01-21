process COSMIC_DOWNLOAD {
    label 'process_medium'
    label 'process_single_thread'

    input:

        val username
        val password
        val download_url_genes
        val download_url_mutants
        val download_url_celllines_genes
        val download_url_celllines_mutants

    output:
    
    path('cosmic_genes.tar')       , emit: cosmic_genes            , optional:true
    path('cosmic_mutations.tar')   , emit: cosmic_mutants          , optional:true
    path('celllines_genes.tar')    , emit: cosmic_celllines_genes  , optional:true
    path('celllines_mutations.tar'), emit: cosmic_celllines_mutants, optional:true
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    def name   = "${username}"

    """

    key=`echo '${username}:${password}' | base64`
    download_key_genes=`curl -H "Authorization: Basic \$key" "${download_url_genes}"`
    curl "\${download_key_genes:8:-2}" --output cosmic_genes.tar

    download_key_mutants=`curl -H "Authorization: Basic \$key" "${download_url_mutants}"`
    curl "\${download_key_mutants:8:-2}" --output cosmic_mutations.tar

    key=`echo '${username}:${password}' | base64`
    download_key_celllines_genes=`curl -H "Authorization: Basic \$key" "${download_url_celllines_genes}"`
    curl "\${download_key_celllines_genes:8:-2}" --output celllines_genes.tar

    download_key_celllines_mutants=`curl -H "Authorization: Basic \$key" "${download_url_celllines_mutants}"`
    curl "\${download_key_celllines_mutants:8:-2}" --output celllines_mutations.tar
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    
    """
    cosmic_genes.tar
    cosmic_mutations.tar
    celllines_genes.tar
    celllines_mutations.tar

    """
}