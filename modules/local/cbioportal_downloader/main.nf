process CBIOPORTAL_DOWNLOAD {
    label 'process_medium'
    label 'process_single_thread'

    input:

        val cbioportal_download_url
        val sample_id

    output:
    
    path('cbioportal_allstudies_data_mutations.txt')      , emit: cbio_mutations
    path('cbioportal_allstudies_data_clinical_sample.txt'), emit: cbio_samples
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    def name   = "${cbioportal_download_url}"

    """

    wget ${cbioportal_download_url}
    cbioportal=`ls *.tar.gz`
    mv \$cbioportal cbioportal.tar.gz 
    tar -xvzf cbioportal.tar.gz
    cat ${sample_id}/data_mutations.txt > cbioportal_allstudies_data_mutations.txt
    cat ${sample_id}/data_clinical_sample.txt | \\
        awk 'BEGIN{FS=OFS="\\t"}{if(\$1!~"#SAMPLE_ID"){gsub("#SAMPLE_ID", "\\nSAMPLE_ID");} print}' | \\
        awk 'BEGIN{FS=OFS="\\t"}{s=0; j=0; \\
        for(i=1;i<=NF;i++){ \\
            if(\$i=="CANCER_TYPE_DETAILED") j=1; if(\$i=="CANCER_TYPE") s=1; \\
        } \\
        if(j==1 && s==0){gsub("CANCER_TYPE_DETAILED", "CANCER_TYPE");} print;}' \\
        > cbioportal_allstudies_data_clinical_sample.txt


    """


    stub:
    def prefix = task.ext.prefix ?: ""
    
    """
    touch cosmic_genes.tsv

    """
}