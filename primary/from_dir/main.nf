

include { PARSE_SEQ_DIR_UNSPRING } from '../../parse_seq_dir_unspring/main.nf'
include { PRIMARY_FROM_READS } from '../from_reads/main.nf'


workflow PRIMARY_FROM_DIR {
    take: 
    inputDir
    dbPathKraken2
    numReads

    main:

    reads = PARSE_SEQ_DIR_UNSPRING(inputDir)

    PRIMARY_FROM_READS(reads, dbPathKraken2, numReads)

    emit: 
    trimmed = PRIMARY_FROM_READS.out.trimmed
    fastqc_trim_html = PRIMARY_FROM_READS.out.fastqc_trim_html
    fastqc_raw_html = PRIMARY_FROM_READS.out.fastqc_raw_html
    multiqc_html = PRIMARY_FROM_READS.out.multiqc_html
}


