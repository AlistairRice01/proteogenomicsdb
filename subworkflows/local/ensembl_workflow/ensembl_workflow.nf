
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//modules for downloading databases from ENSEMBL
include { PYPGATK_ENSEMBL_DOWNLOAD } from '../../../modules/local/pypgatk_ensembl_downloader/main.nf'
include { CAT_CAT as CAT_DNA       } from '../../../modules/nf-core/cat/cat/main.nf'

//modules for optional database generation
include { PYPGATKDNA as PYPGATK_NCRNA       } from '../../../modules/local/pypgatk_dna/main.nf'
include { PYPGATKDNA as PYPGATK_PSEUDOGENES } from '../../../modules/local/pypgatk_dna/main.nf'
include { PYPGATKDNA as PYPGATK_ALRORFS     } from '../../../modules/local/pypgatk_dna/main.nf'

//modules for the generation of the main ENSEMBL database
include { BCFTOOLS_SORT       } from '../../../modules/nf-core/bcftools/sort/main.nf'
include { BCFTOOLS_CONCAT     } from '../../../modules/nf-core/bcftools/concat/main.nf'
include { TABIX_BGZIP         } from '../../../modules/nf-core/tabix/bgzip/main.nf'
include { BCFTOOLS_INDEX      } from '../../../modules/nf-core/bcftools/index/main.nf'
include { CAT_CAT as CAT_VCF  } from '../../../modules/nf-core/cat/cat/main.nf'
include { PYPGATK_VCF         } from '../../../modules/local/pypgatk_vcf/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        RUN ENSEMBLDB WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ENSEMBLDB {

take:

    //inputs from the main workflow 
    ensembl_downloader_config           //channel: contains the config file for both downloading ENSEMBL data
    species_name                        //channel: contains the species name that will be downloaded from ENSEMBL
    ensembl_config                      //channel: contains the config file for generating an ensembl variant peptide database
    altorfs_config                      //channel: contains the config file for generating a peptide database from altorfs data
    pseudogenes_config                  //channel: contains the config file for generating a peptide database from pseudogenes data
    ncrna_config                        //channel: contains the config file for generating a peptide database from ncrna data

main:

    versions_ch     = Channel.empty()
    mixed_databases = Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF ENSEMBL DOWNLOAD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


    //PYPGATK_ENSEMBL uses the species name and downloads files from ENSEMBL
    PYPGATK_ENSEMBL_DOWNLOAD (
        ensembl_downloader_config,
        species_name
    )
    versions_ch = versions_ch.mix(PYPGATK_ENSEMBL_DOWNLOAD.out.versions).collect()

    //creates an empty channel that will then be populated with the cdna files downloaded from ENSEMBL
    cdna_mixed = Channel.empty()
    cdna_mixed = PYPGATK_ENSEMBL_DOWNLOAD.out.cdna.mix(PYPGATK_ENSEMBL_DOWNLOAD.out.ncrna).collect()
        .map { [ [id: 'total_cDNA' ], it ] }

    //CAT_DNA concatinates the cdna filed downloaded from ENSEMBL into a single file
    CAT_DNA (
        cdna_mixed
    )
    versions_ch = versions_ch.mix(CAT_DNA.out.versions_cat).collect()

    //creates an empty channel that will then be populates with the concatenated cdna
    total_cdna = Channel.empty()
    total_cdna = CAT_DNA.out.file_out.collect()
        .map { meta, it ->
            return [it] }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF ENSEMBL DOWNLOAD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF OPTIONAL DATABASE GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.skip_proteome) {

    //creates an empty channel that will then be populated with the protein files downloaded from ENSEMBL
    reference_proteome = Channel.empty()
    reference_proteome = PYPGATK_ENSEMBL_DOWNLOAD.out.protein.collect()

    //adds the protein files downloaded from ENSEMBL to the mixed database channel
    mixed_databases = mixed_databases.mix(reference_proteome).collect()

}

else {
        //bypass the subworkflow
        log.info "proteome skipped."
    }


if (!params.skip_ncrna) {

    //PYPGATK_NCRNA takes the total_cdna channel and the ncrna_config to generate a peptide database
    PYPGATK_NCRNA (
        total_cdna.map { [ [id: 'ncRNA'], it ] },
        ncrna_config
    )
    versions_ch = versions_ch.mix(PYPGATK_NCRNA.out.versions).collect()

    //adds the peptide database generated from ncrna data to the mixed database channel
    mixed_databases = mixed_databases.mix(PYPGATK_NCRNA.out.database).collect()

}

else {
        //bypass the subworkflow
        log.info "ncrna database skipped."
    }


if (!params.skip_pseudogenes) {

    //PYPGATK_PSEUDOGENES takes the total_cdna channel and the pseudogenes_config to generate a peptide database
    PYPGATK_PSEUDOGENES (
        total_cdna.map { [ [id: 'pseudogenes'], it ] },
        pseudogenes_config
    )
    versions_ch = versions_ch.mix(PYPGATK_PSEUDOGENES.out.versions).collect()

    //adds the peptide database generated from pseudogenes data to the mixed database channel
    mixed_databases = mixed_databases.mix(PYPGATK_PSEUDOGENES.out.database).collect()

}

else {
        //bypass the subworkflow
        log.info "pseudogenes database skipped."
    }

    //creates an empty channel that will then be populated with the cdna files downloaded from ENSEMBL
    cdna_database = Channel.empty()
    cdna_database = PYPGATK_ENSEMBL_DOWNLOAD.out.cdna.collect()
        .map { [ [id: 'Altorfs_database'], it ] }

if (!params.skip_altorfs) {

    //PYPGATK_ALTORFS takes the cdna files and the altorfs_config to generate a peptide database
    PYPGATK_ALRORFS (
        cdna_database,
        altorfs_config
    )
    versions_ch = versions_ch.mix(PYPGATK_ALRORFS.out.versions).collect()

    //adds the peptide database generated from altorfs data to the mixed database channel
    mixed_databases = mixed_databases.mix(PYPGATK_ALRORFS.out.database).collect()
}

else {
        //bypass the subworkflow
        log.info "altorfs database skipped."
    }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF OPTIONAL DATABASE GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        START OF MAIN ENSEMBL DATABASE GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
if (!params.skip_ensembl_vcf) {

    //creates an empty channel which will be populated with the vcf file downloaded from ENSEMBL
    ensembl_vcf  = Channel.empty()
    ensembl_vcf = ensembl_vcf.mix(PYPGATK_ENSEMBL_DOWNLOAD.out.vcf).collect()
        .map { [ [id: 'concatenated_vcf' ], it, [] ] }

    //CAT_VCF concatenates the vcf files downloaded from ENSEMBL into a single file 
    CAT_VCF (
        ensembl_vcf
    )
    versions_ch = versions_ch.mix(CAT_VCF.out.versions_cat).collect()

    //creates an empty channel which will then be populated with the concatenated vcf file 
    total_vcf = Channel.empty()
    total_vcf = CAT_VCF.out.file_out.collect()

    //creates an empty channel that will then be populated with the gtf file downloaded from ENSEMBL
    ensembl_gtf = Channel.empty()
    ensembl_gtf = PYPGATK_ENSEMBL_DOWNLOAD.out.gtf.collect()
        .map { [ [id: 'gtf' ], it ] } 

    //PYPGATK takes the concatenated vcf file, the gtf file, and the concatenated cDNA along with the ensembl_config to generate a peptide database
    PYPGATK_VCF (
        total_vcf,
        ensembl_gtf,
        total_cdna.map { [ [id: 'ensembl_vcf'], it ] },
        ensembl_config
    )
    versions_ch = versions_ch.mix(PYPGATK_VCF.out.versions).collect()

    //adds the peptide database generated from the vcf file to the mixed database channel
    mixed_databases = mixed_databases.mix(PYPGATK_VCF.out.database).collect()
}

else {
        //bypass the subworkflow
        log.info "ensembl vcf database skipped."
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        END OF MAIN ENSEMBL DATABASE GENERATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


emit:
    
    //emits to the main workflow
    mixed_databases     //channel: contains the databases generated from this workflow
    versions_ch         //channel: contains versions.yml holding the version information for each of the tools

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/