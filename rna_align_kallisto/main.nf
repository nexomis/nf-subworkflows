// process include
include { KALLISTO_QUANT } from '../../process/kallisto/quant/main.nf'

// workflow include
include { PARSE_SEQ_DIR_UNSPRING } from '../parse_seq_dir_unspring/main.nf'


workflow RNA_ALIGN_KALLISTO {
    take: 
    inputDir
    idx
    main_settings


    main:

    reads = PARSE_SEQ_DIR_UNSPRING(inputDir)

    KALLISTO_QUANT(
      reads,
      idx,
      main_settings
    )
}


