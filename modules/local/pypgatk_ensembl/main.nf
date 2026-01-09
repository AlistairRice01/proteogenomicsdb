process PYPGATK_ENSEMBL {
    label 'process_medium'
    label 'process_single_thread'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f734a98d2a7b9c84d988474c7e71aef2b0abb987b36fb83b17d52bec0d3babe/data' :
        'community.wave.seqera.io/library/pypgatk_pip:456a1305c1d65d3a' }"

    input:
    path ensembl_downloader_config
    val species_taxonomy

    output:

    path "*.pep.all.fa" , emit: protein
    path "*cdna.all.fa" , emit: cdna
    path "*ncrna.fa", emit: ncrna
    path "*.dna*.fa", emit: fasta
    path "*.gtf", emit: gtf
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''

    """
    pypgatk_cli.py ensembl-downloader \\
        --config_file ${ensembl_downloader_config} \\
        --output_directory . \\
        --taxonomy ${species_taxonomy}  \\
        -sv -sc
 
    cat <<-END_VERSIONS > versions.yml
        pypgatk: \$(echo \$(pypgatk --version 2>&1) | sed 's/^pypgatk v//')
    END_VERSIONS
    """

    stub:
  
    """

    cat <<-END_VERSIONS > versions.yml
        pypgatk: \$(echo \$(pypgatk --version 2>&1) | sed 's/^pypgatk v//')
    END_VERSIONS
    """
}