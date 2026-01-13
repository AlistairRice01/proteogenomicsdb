
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//modules involved in alignment of the reads to the reference genome
include { HISAT2_EXTRACTSPLICESITES } from '../../../modules/nf-core/hisat2/extractsplicesites/main.nf'
include { HISAT2_BUILD              } from '../../../modules/nf-core/hisat2/build/main.nf'
include { HISAT2_ALIGN              } from '../../../modules/nf-core/hisat2/align/main.nf'

//modules involved in the generation of a database from a DNA sequence
include { SAMTOOLS_SORT       } from '../../../modules/nf-core/samtools/sort/main.nf'
include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main.nf'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main.nf'
include { SAMTOOLS_FAIDX      } from '../../../modules/nf-core/samtools/faidx/main.nf'
include { GFFCOMPARE          } from '../../../modules/nf-core/gffcompare/main.nf'
include { GFFREAD             } from '../../../modules/nf-core/gffread/main.nf'
include { PYPGATKDNA          } from '../../../modules/local/pypgatk_dna/main.nf'

//modules involved in the generation of a database from a VCF file
include { SAMTOOLS_INDEX       } from '../../../modules/nf-core/samtools/index/main.nf'
include { FREEBAYES            } from '../../../modules/nf-core/freebayes/main.nf'
include { GUNZIP as GUNZIP_VCF } from '../../../modules/nf-core/gunzip/main.nf'
include { PYPGATKCUSTOM        } from '../../../modules/local/pypgatk_custom/main.nf'

//modules involved in reporting
include { MULTIQC } from '../../../modules/nf-core/multiqc/main.nf'
include { FASTQC  } from '../../../modules/nf-core/fastqc/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQDB {

take:

    // inputs from the main workflow
    reads           //channel: contains the samplesheet of rna-seq reads
    annotation      //channel: contains the annotation (gtf) file 
    reference       //channel: contains the reference genome (fasta) file
    config          //channel: contains the vcf database config file
    dna_config      //channel: contains the DNA database config file
    cdna            //channel: contains the cdna reference (fasta) file
    versions     //channel: contains versions.yml holding the version information for each of the tools 

main:

    //creates an empty channel to combine the two databases generated in this workflow
    Channel
        .empty()
        .set { merged_databases }
        

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF THE ALIGNMENT OF READS TO A REFERENCE PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    //applies a map to the gene annotation
    annotation.map { genome ->
        def meta = [ id: 'genome' ]
        return [ meta, genome ] }
    .set { annotation_ch}
    .collect()

    //HISAT2_EXTRACTSPLICESITES is taking the genome annotation and using it to generate a file containing all of the splicesites in the genome
    HISAT2_EXTRACTSPLICESITES (
        annotation_ch
    )
    versions = versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)

    //this is creating an empty channel which will then be populated with the splicesites data
    Channel
        .empty()
        .set { splicesites_ch }  
    splicesites_ch = HISAT2_EXTRACTSPLICESITES.out.txt.collect()

    //applies a map to the genome reference
    reference.map { genome ->
        def meta = [ id: 'genome' ]
        return [ meta, genome ] }
    .set{ reference_ch }
    .collect()

    //HISAT2_BUILD generates hisat2 index files (.ht2) from the fasta reference, genome annotation, and splicesites 
    HISAT2_BUILD (
        reference_ch, 
        annotation_ch, 
        splicesites_ch
    )
    versions = versions.mix(HISAT2_BUILD.out.versions)

    //creates an empty channel which will then be populated with the genome index files
    Channel
        .empty()
        .set { index_ch }
    index_ch = HISAT2_BUILD.out.index.collect()


    //HISAT2 is taking the fastq files and aligning them to a reference genome using the index files and splicesites 
    HISAT2_ALIGN (
        reads, 
        index_ch, 
        splicesites_ch
    )
    versions = versions.mix(HISAT2_ALIGN.out.versions)

    //creates an empty channel which will then be populated with the HISAT2 alignment data (bam files) 
    Channel
        .empty()
        .set { ch_samtools_input }
    ch_samtools_input = HISAT2_ALIGN.out.bam
        .map { meta, bam ->
            def meta1 = [ id: 'bam_1' ]
            return [ meta1, bam ] }
    
    //creates a channel containing the string 'bai' for specifiying the files produced in SAMTOOLS_SORT
    Channel
        .from('bai')
        .set { sort_input }

    //SAMTOOLS_SORT is taking the HISAT2 alignment data and using the genome reference to sort the files
    SAMTOOLS_SORT (
        ch_samtools_input,
        reference_ch,
        sort_input
    )
    versions = versions.mix(SAMTOOLS_SORT.out.versions_samtools)

    //creates an empty channel which will then be populated with the sorted bam files form SAMTOOLS_SORT
    Channel
        .empty()
        .set { ch_bam_sorted }
    ch_bam_sorted = SAMTOOLS_SORT.out.bam


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF THE ALIGNMENT OF READS TO A REFERENCE PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF THE DNA DATABASE GENERATION PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    //STRINGTIE is taking the genome annotation and the HISAT2 Bam file to assemble the RNA-seq alignments into a potential transcript
    STRINGTIE_STRINGTIE (
        ch_bam_sorted,
        annotation
    )
    versions = versions.mix(STRINGTIE_STRINGTIE.out.versions)

    //creates an empty channel that will then be populated with the STRINGTIE outputs (gtf)
    Channel
        .empty()
        .set{ ch_stringtie_gtf }
    ch_stringtie_gtf = STRINGTIE_STRINGTIE.out.transcript_gtf.map { meta, transcript_gtf ->
        [transcript_gtf] }.collect()

    //STRINGTIE_MERGE combines the gtf files produced by StringTie into a single file
    STRINGTIE_MERGE (
        ch_stringtie_gtf, 
        annotation
    )
        versions = versions.mix(STRINGTIE_MERGE.out.versions)

    //creates an empty channel that will then be populated with the output from STRINGTIE_MERGE (gtf)
    Channel
        .empty()
        .set{ ch_stringtie_gtf_merged }
    ch_stringtie_gtf_merged = STRINGTIE_MERGE.out.gtf
        .map { gtf_file ->
            def meta = [ id: 'stringtie_merged' ]
            return [ meta, gtf_file ] }

    //channel containing a string to set SAMTOOLS_FAIDX to true so that it generates a fai file 
    Channel
        .from('true')
        .set { true_ch }

    //SAMTOOLS_FAIDX takes the reference genome (fasta) and generates a fasta index file (fai)
    SAMTOOLS_FAIDX (
        reference_ch,
        [ [], [] ],
        true_ch  
    )
    versions = versions.mix(SAMTOOLS_FAIDX.out.versions_samtools)

    //creates an empty channel that will then be populated with the fasta index (fai) from SAMTOOLS_FAIDX
    Channel
        .empty()
        .set { ch_fasta_fai }
    ch_fasta_fai = SAMTOOLS_FAIDX.out.fai.collect()

    //cerates an empty channel that will combine the genome fasta and genome fasta index channels
    Channel
        .empty()
        .set { ch_fasta_meta_fai }
    ch_fasta_meta_fai = reference_ch.combine(ch_fasta_fai)
        .map { map1, fasta, map2, fai ->
            def meta = [ id: 'fasta', description: 'Genome FASTA with index' ]
            return [ meta, fasta, fai ] }


    //GFFCOMPARE is taking the genome fasta, genome annotation, and comparing it to the output of StringTie to evaluate the assembly
    GFFCOMPARE (
        ch_stringtie_gtf_merged,
        ch_fasta_meta_fai,
        annotation_ch
    )
    versions = versions.mix(GFFCOMPARE.out.versions)

    //creates an empty channel that will then be populated with the gffcompare output (gtf)
    Channel
        .empty()
        .set{ gffcompare_ch }
    gffcompare_ch = GFFCOMPARE.out.annotated_gtf.collect()
            .map { meta, gtf_file ->
            def meta1 = [ id: 'gtf_file', single_end:true ]
            return [ meta1, gtf_file ] }

    //GFFREAD takes the GFFCOMPARE output (gtf) and using the reference genome converts it into transcript sequences
    GFFREAD (
        gffcompare_ch,
        reference
    )
    versions = versions.mix(GFFREAD.out.versions)
    
    //creates an empty channel that will then be populated with the output from GFFREAD (fasta)
    Channel
        .empty()
        .set { gffout }
    gffout = GFFREAD.out.gffread_fasta.collect()

    //PYPGATK takes the GFFREAD output (FASTA) and using the dna_config file generates a database from the fasta transcripts 
    PYPGATKDNA (
        gffout,
        dna_config
    )
    versions = versions.mix(PYPGATKDNA.out.versions)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF THE DNA DATABASE GENERATION PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    //adds the DNA database to the merged_databases_ch
    merged_databases = merged_databases.mix(PYPGATKDNA.out.database).collect()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF THE VCF DATABASE GENERATION PATHWAY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    //SAMTOOLS_INDEX takes the sorted bam file and generates a bam index file (bai)
    SAMTOOLS_INDEX (
        ch_bam_sorted
    )
    versions = versions.mix(SAMTOOLS_INDEX.out.versions)

    //creates an empty channel that will then be populated with the bam index (bai) generated by SAMTOOLS_INDEX
    Channel
        .empty()
        .set { ch_genome_bai }
    ch_genome_bai = SAMTOOLS_INDEX.out.bai
    .map { bam, bai ->
            return [ bai ] }

    //creates an empty channel that will combine the bam and bam index files into a single channel 
    Channel
        .empty()
        .set { ch_genome_bam_bai }
    ch_genome_bam_bai = ch_bam_sorted.combine(ch_genome_bai)
        .map { meta, bam, bai ->
            return [ meta, bam, bai, [], [], [] ] }

    //FREEBAYES takes the bam, bam index, and fasta index files to generate a vcf file 
    FREEBAYES (
    ch_genome_bam_bai,
    reference_ch,
    ch_fasta_fai,
    [ [], [] ],
    [ [], [] ],
    [ [], [] ]
    )
    versions = versions.mix(FREEBAYES.out.versions)

    //creates an empty channel that will be populated with the gzipped vcf file generated by FREEBAYES
    Channel
        .empty()
        .set { freebayes_ch }
    freebayes_ch = FREEBAYES.out.vcf.collect()

    //GUNZIP_VCF unzips the vcf file generated by FREEBAYES
    GUNZIP_VCF (
        freebayes_ch
    )
    versions = versions.mix(GUNZIP_VCF.out.versions)

    //creates an empty channel that will be populated with the unzipped vcf file generated by FREEBAYES
    Channel
        .empty()
        .set { vcf_unzipped }
    vcf_unzipped = GUNZIP_VCF.out.gunzip.collect()

    //applies a map to the cdna file
    cdna.map { genome ->
        def meta = [ id: 'cdna' ]
        return [ meta, genome ] }
    .set{ cdna_ch }
    .collect()

    //PYPGATKCUSTOM takes the vcf file generated by FREEBAYES and using cdna transcripts and the genome annotation it generates a variant peptide database
    PYPGATKCUSTOM (
        vcf_unzipped,
        annotation_ch,
        cdna_ch,
    )
    versions = versions.mix(PYPGATKCUSTOM.out.versions)

/*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF THE VCF DATABASE GENERATION PATHWAY
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

    //takes the database generated with PYPGATK and combines it with the other database
    merged_databases = merged_databases.mix(PYPGATKCUSTOM.out.database).collect()

emit:

    // emits to the main workflow
    merged_databases         //channel: contains the databases generated from this workflow
    versions                 //channel: contains versions.yml holding the version information for each of the tools 

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/