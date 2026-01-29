process COSMIC_DOWNLOAD {
    label 'process_medium'
    label 'process_single_thread'

    input:

        val username
        val password
        val download_url_genes
        val download_url_mutants

    output:
    
    path('cosmic_genes.fasta')     , emit: cosmic_genes            , optional:true
    path('cosmic_mutations.tsv')   , emit: cosmic_mutants          , optional:true
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    def name   = "${username}"

    """
    download_key_genes=`curl -u ${username}:${password} "${download_url_genes}"`
    curl "\${download_key_genes:8:-2}" --output cosmic_genes.tar
    tar -xvf 'cosmic_genes.tar' '*.gz'
    genes_name=`ls *.gz`
    mv \$genes_name cosmic_genes.fasta.gz
    gunzip cosmic_genes.fasta.gz

    download_key_mutants=`curl -u ${username}:${password} "${download_url_mutants}"`
    curl "\${download_key_mutants:8:-2}" --output cosmic_mutations.tar
    tar -xf cosmic_mutations.tar *.gz
    mutations_name=`ls *.gz`
    mv \$mutations_name cosmic_mutations.tsv.gz
    gunzip cosmic_mutations.tsv.gz

    """

    stub:
    def prefix = task.ext.prefix ?: ""
    
    """
    touch cosmic_genes.fasta
    cosmic_mutations.tsv

    """
}