include { KRAKEN2 as KRAKEN2_HOST1; KRAKEN2 as KRAKEN2_HOST2; KRAKEN2 as KRAKEN2_HOST3 } from '../../process/kraken2/main.nf'
include { GZ as GZ1; GZ as GZ2; GZ as GZ3 } from '../../process/gz/main.nf'
include { checkMeta } from '../utils.nf'

workflow ITERATIVE_UNCLASSIFIED_READS_EXTRACTION {

  take:
  inputReads
  inputK2Index

  main:

  def expectedMeta = [
      "class_db_ids": ["String[]", "ArrayList"],
      "class_tool": ["String"]
  ]

  inputReads.map { checkMeta(it, expectedMeta) }
  inputK2Index.map { checkMeta(it) }

  inputReads.map { [it[0].class_db_ids.size() >= 1 ? it[0].class_db_ids[0] : "", it] }
  | filter { it[0] != "" }
  | combine(inputK2Index.map{[it[0].id, it]}, by:0)
  | set {joinInputForK2i1}

  KRAKEN2_HOST1(joinInputForK2i1.map { it[1] }, joinInputForK2i1.map { it[2] })
  KRAKEN2_HOST1.out.unclassified_reads_fastq
  | GZ1
  | concat(inputReads)
  | unique { it[0].id }
  | set {inputReadsFromK1}

  inputReadsFromK1.map { [it[0].class_db_ids.size() >= 2 ? it[0].class_db_ids[1] : "", it] }
  | filter { it[0] != "" }
  | combine(inputK2Index.map{[it[0].id, it]}, by:0)
  | set {joinInputForK2i2}

  KRAKEN2_HOST2(joinInputForK2i2.map { it[1] }, joinInputForK2i2.map { it[2] })
  KRAKEN2_HOST2.out.unclassified_reads_fastq
  | GZ2
  | concat(inputReadsFromK1)
  | unique { it[0].id }
  | set {inputReadsFromK2}

  inputReadsFromK2.map { [it[0].class_db_ids.size() >= 3 ? it[0].class_db_ids[2] : "", it] }
  | filter { it[0] != "" }
  | combine(inputK2Index.map{[it[0].id, it]}, by:0)
  | set {joinInputForK2i3}

  KRAKEN2_HOST3(joinInputForK2i3.map { it[1] }, joinInputForK2i3.map { it[2] })
  KRAKEN2_HOST3.out.unclassified_reads_fastq
  | GZ3
  | concat(inputReadsFromK2)
  | unique { it[0].id }
  | set {inputReadsFromK3}

  emit:
  inputReadsFromK3

}