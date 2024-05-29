
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

process UNSPRING {
  container 'ghcr.io/nexomis/spring:1.1.1'

  label 'cpu_low'
  label 'mem_12G'

  input:
  tuple val(sample_name), path(spring_file, arity: 1)

  output:
  tuple val("${sample_name}"), path("${sample_name}*.fq.gz", arity: 1..2)

  script:
  """
  #!/usr/bin/bash

  spring -d -g -t ${task.cpus} -i ${sample_name}.spring -o ${sample_name}.fq.gz
  if [[ -e ${sample_name}.fq.gz.1 && -e ${sample_name}.fq.gz.2 ]]; then
    # Move (rename) the files
    mv ${sample_name}.fq.gz.1 ${sample_name}_R1.fq.gz
    mv ${sample_name}.fq.gz.2 ${sample_name}_R2.fq.gz
  fi

  """

  stub:
  """
  #!/usr/bin/bash

  touch ${sample_name}.fq.gz

  """
}

workflow PARSE_SEQ_DIR_UNSPRING {
  take:
  inputDir // either (dir [samples]) or just dir

  main:

  reads = PARSE_SEQ_DIR(inputDir)

  reads.spring  
  | UNSPRING
  | set {unspringReads}

  reads.fastq
  | concat(unspringReads)
  | set {allFastq}

  emit:
  allFastq

}
