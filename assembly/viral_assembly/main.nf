// include - process
include { KRAKEN2 as KRAKEN2_HOST1; KRAKEN2 as KRAKEN2_HOST2; KRAKEN2 as KRAKEN2_HOST3 } from '../../../process/kraken2/main.nf'
include { GZ as GZ1; GZ as GZ2; GZ as GZ3 } from '../../../process/gz/main.nf'
include { BOWTIE2_BUILD } from '../../../process/bowtie2/bowtie2-build/main.nf'
include { BOWTIE2 } from '../../../process/bowtie2/mapping/main.nf'
include { SAM_BAM_SORT_IDX } from '../../../process/samtools/convert-sort-index/main.nf'
include { QUAST } from '../../../process/quast/main.nf'
include { SPADES } from '../../../process/spades/main.nf'
include { ABACAS } from '../../../process/abacas/main.nf'
include { EMPTY_FILE } from '../../../process/empty_file/main.nf'

workflow VIRAL_ASSEMBLY {

  take:
  inputs // (id, [meta, reads], k2_index, inputRefgenome)

  main:
  inputReads = inputs.map { [it[0], it[1]] }
  inputK2Index = inputs.filter { it[2] } | map { [it[0], it[2]] }
  inputRefGenome = inputs.filter { it[2] } | map { [it[0], it[3]] }

  inputK2i1 = inputK2Index.map { [it[0], it[1][0]] }
  inputK2i2 = inputK2Index.filter { it[1].size() > 1 } | map { [it[0], it[1][1]] }
  inputK2i3 = inputK2Index.filter { it[1].size() > 2 } | map { [it[0], it[1][2]] }

  joinInputForK2i1 = inputReads.join(inputK2i1, by:0)
  KRAKEN2_HOST1(joinInputForK2i1.map { it[1] }, joinInputForK2i1.map { it[2] })
  KRAKEN2_HOST1.out.unclassified_reads_fastq
  | GZ1
  | map {[it[0].id, it]}
  | concat(inputReads)
  | unique { it[0] }
  | set {inputReadsFromK1}

  joinInputForK2i2 = inputReadsFromK1.join(inputK2i2, by:0)
  KRAKEN2_HOST2(joinInputForK2i2.map { it[1] }, joinInputForK2i2.map { it[2] })
  KRAKEN2_HOST2.out.unclassified_reads_fastq
  | GZ2
  | map {[it[0].id, it]}
  | concat(inputReadsFromK1)
  | unique { it[0] }
  | set {inputReadsFromK2}

  joinInputForK2i3 = inputReadsFromK2.join(inputK2i3, by:0)
  KRAKEN2_HOST3(joinInputForK2i3.map { it[1] }, joinInputForK2i3.map { it[2] })
  KRAKEN2_HOST3.out.unclassified_reads_fastq
  | GZ3
  | map {[it[0].id, it]}
  | concat(inputReadsFromK2)
  | unique { it[0] }
  | set {inputReadsFromK3}

  inputReadsFromK3
  | map {it[1]}
  | branch {
      spades: it[0].assembler == "spades"
      no_assembly: true
    }
  | set { inputForAssembly }

  SPADES(inputForAssembly.spades)

  SPADES.out.scaffolds
  | set { scaffolds }

  scaffolds
  | map {[it[0].id, it]}
  | join(inputRefGenome)
  | set { joinInputForAbacas }

  // TODO SELECT BIGGER SCAFFOLD FOR NO_ABACAS ?? What about multi segment ?

  ABACAS(joinInputForAbacas.map {it[1]}, joinInputForAbacas.map {it[2]})
  ABACAS.out
  | concat(scaffolds)
  | unique { it[0].id }
  | set {finalScaffolds}

  // mapping (with index building and bam sorting)

  finalScaffolds
  | filter { it[0].realign == "yes" }
  | set { toRealignScaffolds }

  BOWTIE2_BUILD(toRealignScaffolds)
  BOWTIE2_BUILD.out.idx
  | map { [it[0].id, it[1]] }
  | set { bwtIdx }

  joinInputforBt2 = inputReadsFromK3.join(bwtIdx)

  BOWTIE2(joinInputforBt2.map {it[1]}, joinInputforBt2.map {it[2]})

  rawSAM=BOWTIE2.out.sam
  SAM_BAM_SORT_IDX(rawSAM)
  SAM_BAM_SORT_IDX.out.bam
  | toSortedList {it[0].id}
  | flatMap
  | set { sortedBAM }

  emptyFile = EMPTY_FILE()

  sortedBAM
  | map { [it[0].id, it[1], it[2]] }
  | set { bamForQuast }

  finalScaffolds
  | map { [it[0].id, it[0], it[1]] }
  | join(inputRefGenome, by: 0, remainder: true)
  | combine(emptyFile)
  | map {
    if (it[3]) {
      return [it[0], it[1], it[2], it[3]]
    } else {
      return [it[0], it[1], it[2], it[4]]
    }
  }
  | join(bamForQuast, by: 0, remainder: true)
  | combine(emptyFile)
  | map {
      if (it[4]) {
        return [it[1], it[2], it[3], it[4], it[5]]
      } else {
        return [it[1], it[2], it[3], it[5], it[5]]
      }
    }
| set {inputforQuast}
  
  QUAST(inputforQuast)

  emit:
  scaffolds = finalScaffolds
  alignments = sortedBAM

}
