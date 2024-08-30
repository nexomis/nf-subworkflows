// include - process
include { KRAKEN2 as KRAKEN2_HOST } from '../../../process/kraken2/main.nf'
include { BOWTIE2_BUILD } from '../../../process/bowtie2/bowtie2-build/main.nf'
include { BOWTIE2 } from '../../../process/bowtie2/mapping/main.nf'
include { SAM_BAM_SORT_IDX } from '../../../process/samtools/convert-sort-index/main.nf'
include { QUAST } from '../../../process/quast/main.nf'
include { SPADES } from '../../../process/spades/main.nf'
include { ABACAS } from '../../../process/abacas/main.nf'


workflow VIRAL_ASSEMBLY {
  // TODO: if usefull manage compression (.gz) of intermediate and final fasta ?

  take:
  inputReads
  inputK2Index
  inputRefGenome

  main:

  inputReads
  | toSortedList {it[0].id}
  | flatMap
  | set {sortedReads}

  //sortedReads = Channel.fromList(sortedReadsList)
  
  sortedReads
  | branch {
      filter: it[0].host_removal == "yes"
      no_filter: true
  }
  | set { inputForK2 }

  
  inputK2Index
  | filter { it[0].host_removal == "yes" }
  | toSortedList {it[0].id}
  | flatMap
  | map { it[1] }
  | set { sortedInputK2Index }
    
  KRAKEN2_HOST(inputForK2.filter, sortedInputK2Index)
  KRAKEN2_HOST.out.unclassified_reads_fastq
  | concat(inputForK2.no_filter)
  | toSortedList {it[0].id}
  | flatMap
  | set {filteredInputReads}

  filteredInputReads
  | branch {
      spades: it[0].assembler == "spades"
      no_assembly: true
    }
  | set { inputForAssembly }

  SPADES(inputForAssembly.spades)
  SPADES.out.scaffolds
  | toSortedList {it[0].id}
  | flatMap
  | set { scaffolds }

  scaffolds
  | branch {
      abacas: it[0].abacas == "yes"
      no_abacas: true
    }
  | set { scaffoldsForAbacas }

  inputRefGenome
  | filter { it[0].abacas == "yes" }
  | toSortedList {it[0].id}
  | flatMap
  | map {it[1]}
  | set { sortedRefGenome }

  // TODO SELECT BIGGER SCAFFOLD FOR NO_ABACAS ?? What about multi segment ?

  ABACAS(scaffoldsForAbacas.abacas, sortedRefGenome)
  ABACAS.out.scaffolds
  | concat(scaffoldsForAbacas.no_abacas)
  | toSortedList {it[0].id}
  | flatMap
  | set {finalScaffolds}


  // mapping (with index building and bam sorting)

  finalScaffolds
  | filter { it[0].realign == "yes" }
  | set { toRealignScaffolds }

  filteredInputReads
  | filter { it[0].realign == "yes" }
  | set { readsForBowtie2 }

  //toRealignScaffolds
  //| view

  BOWTIE2_BUILD(toRealignScaffolds)
  BOWTIE2_BUILD.out.idx
  | toSortedList {it[0].id}
  | flatMap
  | map { it[1] }
  | set { bwtIdx }

  bwtIdx
  | view

  readsForBowtie2
  | view

  BOWTIE2(readsForBowtie2, bwtIdx)

  rawSAM=BOWTIE2.out.sam
  SAM_BAM_SORT_IDX(rawSAM)
  SAM_BAM_SORT_IDX.out.bam
  | toSortedList {it[0].id}
  | flatMap
  | set { sortedBAM }
  
  // QUAST(scaffoldsAssembled, sortedBAM, inputRefGenome)
  // QUAST(scaffoldsAssembled, sortedBAM, inputRefGenome)
  // QUAST(scaffoldsAssembled, null, null)

  emit:
  scaffolds = finalScaffolds
  alignments = sortedBAM
  
}
