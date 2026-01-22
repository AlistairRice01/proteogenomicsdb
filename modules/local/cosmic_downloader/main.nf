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
    
    path('cosmic_genes.tsv')       , emit: cosmic_genes            , optional:true
    path('cosmic_mutations.tsv')   , emit: cosmic_mutants          , optional:true
    path('celllines_genes.tsv')    , emit: cosmic_celllines_genes  , optional:true
    path('celllines_mutations.tsv'), emit: cosmic_celllines_mutants, optional:true
    
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
    tar -xvf 'cosmic_genes.tar' '*.gz'
    genes_name=`ls *.gz`
    mv \$genes_name cosmic_genes.tsv.gz
    gunzip cosmic_genes.tsv.gz

        download_key_mutants=`curl -H "Authorization: Basic \$key" "${download_url_mutants}"`
    curl "\${download_key_mutants:8:-2}" --output cosmic_mutations.tar
    tar -xf cosmic_mutations.tar *.gz
    mutations_name=`ls *.gz`
    mv \$mutations_name cosmic_mutations.tsv.gz
    gunzip cosmic_mutations.tsv.gz

    key=`echo '${username}:${password}' | base64`
    download_key_celllines_genes=`curl -H "Authorization: Basic \$key" "${download_url_celllines_genes}"`
    curl "\${download_key_celllines_genes:8:-2}" --output celllines_genes.tar
    tar -xf celllines_genes.tar *.gz
    celllines_genes_name=`ls *.gz`
    mv \$celllines_genes_name celllines_genes.tsv.gz
    gunzip celllines_genes.tsv.gz

    download_key_celllines_mutants=`curl -H "Authorization: Basic \$key" "${download_url_celllines_mutants}"`
    curl "\${download_key_celllines_mutants:8:-2}" --output celllines_mutations.tar
    file=`tar -xf celllines_mutations.tar *.gz`
    celllines_mutations_name=`ls *.gz`
    mv \$celllines_mutations_name celllines_mutations.tsv.gz
    gunzip celllines_mutations.tsv.gz

    """


    stub:
    def prefix = task.ext.prefix ?: ""
    
    """
    touch cosmic_genes.tsv
    touch cosmic_mutations.tsv
    touch celllines_genes.tsv
    touch celllines_mutations.tsv

    """
}