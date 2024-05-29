// process include
include { FASTP } from '../../process/fastp/main.nf'

// workflow include
include { PARSE_SEQ_DIR_UNSPRING } from '../parse_seq_dir_unspring/main.nf'

workflow PRIMARY {
    take: 
    inputDir
    fastpArgs

    main:

    reads = PARSE_SEQ_DIR_UNSPRING(inputDir)

    trimmedReads = FASTP(
      reads,
      fastpArgs
    )

    emit: 
    trimmed = trimmedReads

}
