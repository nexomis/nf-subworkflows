// process include
include { FASTP } from '../../process/fastp/main.nf'
include { KRAKEN2 } from '../../process/kraken2/main.nf'

include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIMMED } from '../../process/fastqc/main.nf'

process PRIMARY_MULTIQC {
  container 'multiqc/multiqc:v1.21' // version above are bugged

  label 'cpu_x1'
  label 'mem_8G'

  input:
  path raw_zips, stageAs: 'raw/*', arity: '1..*'
  path trimmed_zips, stageAs: 'trimmed/*', arity: '1..*'
  path kraken2_reports, arity: '1..*'
  path conf_yml, arity: 1, stageAs: "multiqc_config.yaml"

  output:
  path("primary_multiqc.html", arity: 1)

  script:
  """
  #!/usr/bin/bash

  multiqc -c $conf_yml --no-data-dir -n primary_multiqc.html .

  """

  stub:
  """
  #!/usr/bin/bash

  touch primary_multiqc.html

  """

}

// workflow include
include { PARSE_SEQ_DIR_UNSPRING } from '../parse_seq_dir_unspring/main.nf'

workflow PRIMARY {
    take: 
    inputDir
    dbPathKraken2

    main:

    reads = PARSE_SEQ_DIR_UNSPRING(inputDir)

    trimmedReads = FASTP(
      reads
    )

    FASTQC_RAW(reads)

    FASTQC_TRIMMED(trimmedReads)

    Channel.fromPath(moduleDir + "/multiqc.yml")
    | set {multiqcYml}

    FASTQC_RAW.out.zip
    | map {it[1]}
    | collect
    | set {fastqcRaw}
    
    FASTQC_TRIMMED.out.zip
    | map {it[1]}
    | collect
    | set {fastqcTrimmed}

    view(trimmedReads)

    KRAKEN2(
      trimmedReads,
      dbPathKraken2
    )

    KRAKEN2.out.report
    | map {it[1]}
    | collect
    | set {kraken2Reports}

    PRIMARY_MULTIQC(
      fastqcRaw,
      fastqcTrimmed,
      kraken2Reports,
      multiqcYml
    )

    emit: 
    trimmed = trimmedReads

}


