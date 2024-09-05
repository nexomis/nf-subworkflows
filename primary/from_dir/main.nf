

include { PARSE_SEQ_DIR_UNSPRING } from '../../parse_seq_dir_unspring/main.nf'
include { PRIMARY_FROM_READS } from '../from_reads/main.nf'


workflow PRIMARY_FROM_DIR {
    take: 
    inputDir
    dbPathKraken2

    main:

    reads = PARSE_SEQ_DIR_UNSPRING(inputDir)

    PRIMARY_FROM_READS(reads, dbPathKraken2)

    emit: 
    trimmed = PRIMARY_FROM_READS.out.trimmed

}


