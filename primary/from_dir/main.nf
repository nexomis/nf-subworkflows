

include { PARSE_SEQ_DIR_UNSPRING } from '../../parse_seq_dir_unspring/main.nf'
include { PRIMARY } from '../main.nf'


workflow PRIMARY_FROM_DIR {
    take: 
    inputDir
    dbPathKraken2
    numReads

    main:

    reads = PARSE_SEQ_DIR_UNSPRING(inputDir)

    PRIMARY(reads, dbPathKraken2, numReads)

    emit: 
    trimmed = PRIMARY.out.trimmed
    fastqc_trim_html = PRIMARY.out.fastqc_trim_html
    fastqc_raw_html = PRIMARY.out.fastqc_raw_html
    multiqc_html = PRIMARY.out.multiqc_html
    kraken2_report = PRIMARY.out.kraken2_report
    kraken2_output = PRIMARY.out.kraken2_output
    class_report = PRIMARY.out.class_report
}


