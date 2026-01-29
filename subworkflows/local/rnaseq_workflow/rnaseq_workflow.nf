
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC_WORKFLOW } from '../fastqc_workflow/fastqc_workflow.nf'

//modules involved in alignment of the reads to the reference genome
include { HISAT2_EXTRACTSPLICESITES } from '../../../modules/nf-core/hisat2/extractsplicesites/main.nf'
include { HISAT2_BUILD              } from '../../../modules/nf-core/hisat2/build/main.nf'
include { HISAT2_ALIGN              } from '../../../modules/nf-core/hisat2/align/main.nf'

//modules involved in the generation of a database from a DNA sequence
include { SAMTOOLS_SORT       } from '../../../modules/nf-core/samtools/sort/main.nf'
include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main.nf'
include { SAMTOOLS_FAIDX      } from '../../../modules/nf-core/samtools/faidx/main.nf'
include { GFFCOMPARE          } from '../../../modules/nf-core/gffcompare/main.nf'
include { GFFREAD             } from '../../../modules/nf-core/gffread/main.nf'
include { PYPGATKDNA          } from '../../../modules/local/pypgatk/dnaseq_to_proteindb/main.nf'

//modules involved in the generation of a database from a VCF file
include { FREEBAYES            } from '../../../modules/nf-core/freebayes/main.nf'
include { GUNZIP as GUNZIP_VCF } from '../../../modules/nf-core/gunzip/main.nf'
include { PYPGATK_VCF         } from '../../../modules/local/pypgatk/vcf_to_proteindb/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQDB {

take:

    //inputs from the main workflow
    samplesheet            //channel: [val(meta), [ samplesheet ] ]
    annotation             //channel: /path/to/genome annotation
    reference              //channel: /path/to/reference genome
    config                 //channel: /path/to/custom config file
    dna_config             //channel: /path/to/dna config file
    cdna                   //channel: /path/to/cdna transcripts
    faidx_get_genome_sizes //boolean: whether FAIDX gets genome sizes
    samtools_sort_index    //string: index sorting type     
    multiqc_config         //channel: /path/to/multiqc config


main:

    //creates empty channels for tool versions, the peptide database, and the multiqc report
    merged_databases_ch = Channel.empty()
    multiqc_report_ch   = Channel.empty()
    versions_ch         = Channel.empty()

    //FASTQC_WORKFLOW takes the sample sheet and mutltiqc and runs the reads through fastqc, trimgalore, and multiqc
    FASTQC_WORKFLOW (
    samplesheet,
    multiqc_config
    )
    //extract tool versions and the multiqc report from FASTQC_WOKFLOW
    versions_ch       = versions_ch.mix(FASTQC_WORKFLOW.out.versions_ch).collect()
    multiqc_report_ch = FASTQC_WORKFLOW.out.multiqc_report_ch.collect()
        
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF THE ALIGNMENT OF READS TO A REFERENCE PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    //creating empty channels used in the alignment pathway
    splicesites_ch   = Channel.empty()
    index_ch         = Channel.empty()
    bam_alignment_ch = Channel.empty()
    bam_sorted_ch    = Channel.empty()

    //applies a map to the genome annotation
    annotation.map { genome ->
        def meta = [ id: 'Genome' ]
        return [ meta, genome ] }
    .set { annotation_ch}
    .collect()

    //HISAT2_EXTRACTSPLICESITES takes the genome annotation and using it to 
    //generate a file containing all of the splicesites in the genome
    HISAT2_EXTRACTSPLICESITES (
        annotation_ch
    )
    versions_ch = versions_ch.mix(HISAT2_EXTRACTSPLICESITES.out.versions) 
    splicesites_ch = HISAT2_EXTRACTSPLICESITES.out.txt.collect()

    //applies a map to the genome reference
    reference.map { genome ->
        def meta = [ id: 'genome' ]
        return [ meta, genome ] }
    .set{ reference_ch }
    .collect()

    //HISAT2_BUILD generates hisat2 index files (.ht2) 
    //from the fasta reference, genome annotation, and splicesites 
    HISAT2_BUILD (
        reference_ch, 
        annotation_ch, 
        splicesites_ch
    )
    versions_ch = versions_ch.mix(HISAT2_BUILD.out.versions)
    index_ch = HISAT2_BUILD.out.index.collect()


    //HISAT2 is taking the fastq files and aligning them to a 
    //reference genome using the index files and splicesites 
    HISAT2_ALIGN (
        samplesheet, 
        index_ch, 
        splicesites_ch
    )
    versions_ch = versions_ch.mix(HISAT2_ALIGN.out.versions)
    bam_alignment_ch = HISAT2_ALIGN.out.bam
        .map { meta, bam ->
            def meta1 = [ id: 'Bam_1' ]
            return [ meta1, bam ] }

    
    //SAMTOOLS_SORT is taking the HISAT2 alignment 
    //data and using the genome reference to sort the files
    SAMTOOLS_SORT (
        bam_alignment_ch,
        reference_ch,
        samtools_sort_index
    )
    versions_ch = versions_ch.mix(SAMTOOLS_SORT.out.versions_samtools)
    bam_sorted_ch = SAMTOOLS_SORT.out.bam



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF THE ALIGNMENT OF READS TO A REFERENCE PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF THE DNA DATABASE GENERATION PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//conditional execution based on if skip_dnaseq is true or false
if (!params.skip_dnaseq) {

    //creating empty channels used in the dna database generation pathway
    stringtie_gtf_ch        = Channel.empty()
    fasta_fai_ch            = Channel.empty()
    gffcompare_ch           = Channel.empty()
    gffout                  = Channel.empty()

    //STRINGTIE is taking the genome annotation and the HISAT2 Bam file 
    //to assemble the RNA-seq alignments into a potential transcript
    STRINGTIE_STRINGTIE (
        bam_sorted_ch,
        annotation
    )
    versions_ch = versions_ch.mix(STRINGTIE_STRINGTIE.out.versions)
    stringtie_gtf_ch = STRINGTIE_STRINGTIE.out.transcript_gtf

    //SAMTOOLS_FAIDX takes the reference genome (fasta) 
    //and generates a fasta index file (fai)
    SAMTOOLS_FAIDX (
        reference_ch,
        [ [], [] ],
        faidx_get_genome_sizes  
    )
    versions_ch = versions_ch.mix(SAMTOOLS_FAIDX.out.versions_samtools)
    fasta_fai_ch = SAMTOOLS_FAIDX.out.fai.collect()

    //cerates an empty channel that will combine the 
    //genome fasta and genome fasta index channels
    fasta_meta_fai_ch = Channel.empty()
    fasta_meta_fai_ch = reference_ch.combine(fasta_fai_ch)
        .map { map1, fasta, map2, fai ->
            def meta = [ id: 'fasta', description: 'Genome FASTA with index' ]
            return [ meta, fasta, fai ] }


    //GFFCOMPARE is taking the genome fasta, genome annotation, and the 
    //output of StringTie to evaluate the assembly
    GFFCOMPARE (
        stringtie_gtf_ch,
        fasta_meta_fai_ch,
        annotation_ch
    )
    versions_ch = versions_ch.mix(GFFCOMPARE.out.versions)
    gffcompare_ch = GFFCOMPARE.out.annotated_gtf.collect()
            .map { meta, gtf_file ->
            def meta1 = [ id: 'gtf_file', single_end:true ]
            return [ meta1, gtf_file ] }

    //GFFREAD takes the GFFCOMPARE output (gtf) and using 
    //the reference genome converts it into transcript sequences
    GFFREAD (
        gffcompare_ch,
        reference
    )
    versions_ch = versions_ch.mix(GFFREAD.out.versions)
    gffout = GFFREAD.out.gffread_fasta.collect()

    //PYPGATK takes the GFFREAD output (FASTA) and using the dna_config 
    //file generates a database from the fasta transcripts 
    PYPGATKDNA (
        gffout,
        dna_config
    )
    versions_ch = versions_ch.mix(PYPGATKDNA.out.versions)

    //adds the DNA database to the merged_databases_ch
    merged_databases_ch = merged_databases_ch.mix(PYPGATKDNA.out.database).collect()

}

else {
        //bypass the subworkflow
        log.info "rna-seq dna database skipped."
    }
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF THE DNA DATABASE GENERATION PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF THE VCF DATABASE GENERATION PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//conditional execution based on if skip_vcf is true or false
if (!params.skip_vcf) {

        
    genome_bai_ch     = Channel.empty()
    genome_bam_bai_ch = Channel.empty()
    freebayes_ch      = Channel.empty()

    //take the index (bai) file generated by SAMTOOLS_SORT
    genome_bai_ch = SAMTOOLS_SORT.out.bai
        .map { meta, bai ->
            return [ bai ] }

    //combines the bam and bam index files into a single channel
    genome_bam_bai_ch = bam_sorted_ch.combine(genome_bai_ch)
        .map { meta, bam, bai ->
            return [ meta, bam, bai, [], [], [] ] }

    //FREEBAYES takes the bam, bam index, and fasta index files to generate a vcf file 
    FREEBAYES (
    genome_bam_bai_ch,
    reference_ch,
    fasta_fai_ch,
    [ [], [] ],
    [ [], [] ],
    [ [], [] ]
    )
    versions_ch = versions_ch.mix(FREEBAYES.out.versions)
    freebayes_ch = FREEBAYES.out.vcf.collect()

    //GUNZIP_VCF unzips the vcf file generated by FREEBAYES
    GUNZIP_VCF (
        freebayes_ch
    )
    versions_ch = versions_ch.mix(GUNZIP_VCF.out.versions)

    //creates an empty channel that will be populated with the unzipped vcf file generated by FREEBAYES
    vcf_unzipped = Channel.empty()
    vcf_unzipped = GUNZIP_VCF.out.gunzip.collect()

    //applies a map to the cdna file
    cdna.map { genome ->
        def meta = [ id: 'cdna' ]
        return [ meta, genome ] }
    .set{ cdna_ch }
    .collect()

    //PYPGATKCUSTOM takes the vcf file generated by FREEBAYES and using cdna 
    //transcripts and the genome annotation it generates a variant peptide database
    PYPGATK_VCF (
        vcf_unzipped,
        annotation_ch,
        cdna_ch,
        config
    )
    versions_ch = versions_ch.mix(PYPGATK_VCF.out.versions)

    //takes the database generated with PYPGATK and combines it with the other database
    merged_databases_ch = merged_databases_ch.mix(PYPGATK_VCF.out.database).collect()

}

else {
    //bypass the subworkflow
    log.info "rna-seq vcf database skipped."
}

/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF THE VCF DATABASE GENERATION PATHWAY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


emit:

    //emits to the main workflow
    merged_databases_ch //channel: [ val(meta), [ database ] ]
    versions_ch         //channel: [ path(versions.yml) ]
    multiqc_report_ch   //channel: /path/to/multiqc_report.html

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/