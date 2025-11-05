include { SPRING_DECOMPRESS } from '../../process/spring/decompress/main.nf'
include { checkMeta } from '../utils.nf'

workflow UNSPRING_READS {

  take:
  inputReads

  main:

  // Validate metadata
  inputReads.map { it -> checkMeta(it) }

  inputReads
  | branch { it -> 
    spring: it[0].read_type == "spring"
    fastq: true
  }
  | set {reads}

  reads.spring
  | SPRING_DECOMPRESS
  | map{ it -> 
    if (it[1].size() == 1) {
        it[0].read_type = "SR"
    } else {
        it[0].read_type = "PE"
    }
    return(it)
  }
  | set {unspringReads}

  reads.fastq
  | concat(unspringReads)
  | set {allFastq}

  emit:
  allFastq

}
