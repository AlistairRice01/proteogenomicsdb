process PYPGATK_VCF {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_single_thread'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f734a98d2a7b9c84d988474c7e71aef2b0abb987b36fb83b17d52bec0d3babe/data' :
        'community.wave.seqera.io/library/pypgatk_pip:456a1305c1d65d3a' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta), path(gtf)
    tuple val(meta), path(dna)
    file pypgatk_config

    output:
    tuple val(meta), path("*.fa"), emit: database
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def name = "${prefix}.fa"

    """
    pypgatk_cli.py vcf-to-proteindb \\
        --config_file ${pypgatk_config} \\
        --vcf ${vcf} \\
        --input_fasta ${dna} \\
        --gene_annotations_gtf ${gtf} \\
        --output_proteindb ${name} \\
 
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
    /*   
pypgatk_cli.py vcf-to-proteindb \\
        --config_file ${pypgatk_config} \\
        --vcf ${vcf.baseName}_changedChrNames.vcf \\
        --input_fasta ${dna} \\
        --gene_annotations_gtf ${gtf} \\
        --output_proteindb ${name} \\
        --translation_table 1 \\
        --mito_translation_table 2 \\
        --protein_prefix 'testvar_' \\
        --annotation_field_name '' \\
        --af_field '' \\
        --af_threshold 0.01 \\
        --transcript_str 'FEATURE' \\
        --consequence_str 'CONSEQUENCE' \\
        --accepted_filters 'PASS' \\
        --exclude_consequences 'downstream_gene_variant, upstream_gene_variant, intergenic_variant, intron_variant, synonymous_variant, regulatory_region_variant' \\
        --include_consequences 'all' \\
        --biotype_str transcript_biotype \\
        --exclude_biotypes '' \\
        --include_biotypes 'protein_coding,polymorphic_pseudogene,non_stop_decay,nonsense_mediated_decay,IG_C_gene,IG_D_gene,IG_J_gene,IG_V_gene,TR_C_gene,TR_D_gene,TR_J_gene,TR_V_gene,TEC,mRNA' \\

 python pypgatk_cli.py vcf-to-proteindb
   --vcf sample.vcf
   --input_fasta transcripts.fa
   --gene_annotations_gtf genes.gtf
   --annotation_field_name ''
   --output_proteindb var_peptides.fa

        --vcf ${vcf.baseName}_changedChrNames.vcf \\


    awk 'BEGIN{FS=OFS="\t"}{if(\$1=="chrM") \$1="MT"; gsub("chr","",\$1); print}' \\
        ${vcf} > ${vcf.baseName}_changedChrNames.vcf

*/