process PYPGATK_CBIOPORTAL_DOWNLOAD {
    label 'process_medium'
    label 'process_single_thread'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f734a98d2a7b9c84d988474c7e71aef2b0abb987b36fb83b17d52bec0d3babe/data' :
        'community.wave.seqera.io/library/pypgatk_pip:456a1305c1d65d3a' }"

    input:
        tuple val(meta), path(cbioportal_config)
        val cbioportal_study_id
    
    output:
    path('cbioportal_allstudies_data_mutations.txt'), emit: cbio_mutations
    path('cbioportal_allstudies_data_clinical_sample.txt'), emit: cbio_samples
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def name = "${prefix}"

    if (cbioportal_study_id == "all")
        """
        git clone https://github.com/cBioPortal/datahub.git .
        git lfs install --local --skip-smudge
        git lfs pull -I public --include "data*clinical*sample.txt"
        git lfs pull -I public --include "data_mutations_mskcc.txt"
        cat public/*/data_mutations_mskcc.txt > cbioportal_allstudies_data_mutations.txt
        cat public/*/*data*clinical*sample.txt | \\
            awk 'BEGIN{FS=OFS="\\t"}{if(\$1!~"#SAMPLE_ID"){gsub("#SAMPLE_ID", "\\nSAMPLE_ID");} print}' | \\
            awk 'BEGIN{FS=OFS="\\t"}{s=0; j=0; \\
                for(i=1;i<=NF;i++){ \\
                    if(\$i=="CANCER_TYPE_DETAILED") j=1; \\
                    if(\$i=="CANCER_TYPE") s=1; \\
                } \\
                if(j==1 && s==0){ \\
                    gsub("CANCER_TYPE_DETAILED", "CANCER_TYPE"); \\
                } \\
                print; \\
            }' \\
            > cbioportal_allstudies_data_clinical_sample.txt \\

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pypgatk: \$(echo \$(pypgatk --version 2>&1) | sed 's/^pypgatk v//')
        END_VERSIONS
        """
    else
        """
        pypgatk_cli.py cbioportal-downloader \\
            --config_file ${cbioportal_config} \\
            --download_study ${cbioportal_study_id} \\
            --output_directory ${name} \\
        
        tar -xvf ${name}/${cbioportal_study_id}.tar
        cat ${cbioportal_study_id}/data_mutations.txt > cbioportal_allstudies_data_mutations.txt
        cat ${cbioportal_study_id}/data_clinical_sample.txt | \\
            awk 'BEGIN{FS=OFS="\\t"}{if(\$1!~"#SAMPLE_ID"){gsub("#SAMPLE_ID", "\\nSAMPLE_ID");} print}' | \\
            awk 'BEGIN{FS=OFS="\\t"}{s=0; j=0; \\
            for(i=1;i<=NF;i++){ \\
                if(\$i=="CANCER_TYPE_DETAILED") j=1; if(\$i=="CANCER_TYPE") s=1; \\
            } \\
            if(j==1 && s==0){gsub("CANCER_TYPE_DETAILED", "CANCER_TYPE");} print;}' \\
            > cbioportal_allstudies_data_clinical_sample.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            pypgatk: \$(echo \$(pypgatk --version 2>&1) | sed 's/^pypgatk v//')
        END_VERSIONS
        """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pypgatk: \$(echo \$(pypgatk --version 2>&1) | sed 's/^pypgatk v//')
    END_VERSIONS
    """
}
