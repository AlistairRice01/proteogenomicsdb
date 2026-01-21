process COSMIC_DOWNLOAD {
    label 'process_medium'
    label 'process_single_thread'

    input:
        val download_url
        val username
        val password
        val output_file
    
    output:
    path('*.*')         , emit: cosmic_file

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    def name   = "${download_url}"
    
    """
    key=`echo '${username}:${password}' | base64`
    download_key=`curl -H "Authorization: Basic \$key" "${download_url}"`
    curl "\${download_key:8:-2}" --output ${output_file}

    """

    stub:
    def prefix = task.ext.prefix ?: ""
    
    """
    touch cosmic_file.txt

    """
}
