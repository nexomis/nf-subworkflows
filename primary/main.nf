// process include
include { FASTP } from '../../process/fastp/main.nf'
include { KRAKEN2 } from '../../process/kraken2/main.nf'
include { RECENTRIFUGE } from '../../process/recentrifuge/main.nf'
include { SEQTK_SAMPLE as SEQTK_SAMPLE_RAW; SEQTK_SAMPLE as SEQTK_SAMPLE_TRIMMED } from '../../process/seqtk/sample/main.nf'
include { SLIMFASTQ_DECOMPRESS } from '../../process/slimfastq/decompress/main.nf'
include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIMMED } from '../../process/fastqc/main.nf'

process PRIMARY_MULTIQC {
  container 'multiqc/multiqc:v1.21' // version above are bugged

  label 'cpu_x1'
  label 'mem_8G'

  input:
  path raw_zips, stageAs: 'raw/*', arity: '1..*'
  path trimmed_zips, stageAs: 'trimmed/*', arity: '1..*'
  path kraken2_reports, arity: '1..*'
  path fastp_reports, arity: '1..*'
  path conf_yml, arity: 1, stageAs: "multiqc_config.yaml"

  output:
  path("primary_multiqc.html", arity: 1)

  script:
  """
  #!/usr/bin/bash

  multiqc -c $conf_yml --exclude general_stats --no-data-dir -n primary_multiqc.html .

  """

  stub:
  """
  #!/usr/bin/bash

  touch primary_multiqc.html

  """

}

workflow PRIMARY {
  take: 
  inputReads
  dbPathKraken2
  taxDir // got with retaxdump from recentrifuge
  numReads

  main:

  inputReads
  | branch {
    sfq: it[0].read_type == "sfq"
    fastq: it[0].read_type != "sfq"
  }
  | set { inputs }

  SLIMFASTQ_DECOMPRESS(inputs.sfq)
  | map {
    if (it[1].size()==1) {
      it[0].read_type = "SR"
    } else {
      it[0].read_type = "PE"
    }
    return it
  }
  | concat(inputs.fastq)
  | set { inputFastq }


  FASTP(inputFastq)

  trimmedReads = FASTP.out.reads

  SEQTK_SAMPLE_RAW(inputFastq, numReads)
  | set {subInputFastq}
  SEQTK_SAMPLE_TRIMMED(trimmedReads, numReads)
  | set {subTrimmedReads}

  FASTQC_RAW(subInputFastq)

  FASTQC_TRIMMED(subTrimmedReads)

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

  FASTP.out.report  
  | map {it[1]}
  | collect
  | set {fastpReports}

  KRAKEN2(
    subTrimmedReads,
    dbPathKraken2
  )

  KRAKEN2.out.report
  | map {it[1]}
  | collect
  | set {kraken2Reports}

  KRAKEN2.out.output
  | map {it[1]}
  | collect
  | set {kraken2Outputs}

  RECENTRIFUGE(kraken2Outputs, taxDir)

  PRIMARY_MULTIQC(
    fastqcRaw,
    fastqcTrimmed,
    kraken2Reports,
    fastpReports,
    multiqcYml
  )

  emit:
  trimmed = trimmedReads
  fastqc_trim_html = FASTQC_TRIMMED.out.html
  fastqc_raw_html = FASTQC_RAW.out.html
  multiqc_html = PRIMARY_MULTIQC.out
  kraken2_report = KRAKEN2.out.report
  kraken2_output = KRAKEN2.out.output
  class_report = RECENTRIFUGE.out.html

}
