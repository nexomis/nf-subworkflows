
/*

Parse sequence directory and uncompress file if necessary.

The Inputdir can be either 
(directory)
OR
(directory, sample names) in that case only reads matching sample names are
  treated

A subworkflow without process.
It takes directory as input

It gives list of sequences files with sample name as output:
paired: [sample_name [file R1, file R2]]
single: [sample_name [file R1]]

*/

include { PARSE_SEQ_DIR } from '../parse_seq_dir/main.nf'
include { SPRING_DECOMPRESS } from '../../process/spring/decompress/main.nf'

workflow PARSE_SEQ_DIR_UNSPRING {
  take:
  inputDir

  main:

  reads = PARSE_SEQ_DIR(inputDir)

  reads.spring  
  | SPRING_DECOMPRESS
  | set {unspringReads}

  reads.fastq
  | concat(unspringReads)
  | set {allFastq}

  emit:
  allFastq

}
