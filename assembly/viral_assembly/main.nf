// include - process
include { KRAKEN2 as KRAKEN2_HOST } from '../../../process/kraken2/main.nf'
include { BOWTIE2_BUILD } from '../../../process/bowtie2/bowtie2-build/main.nf'
include { BOWTIE2 } from '../../../process/bowtie2/mapping/main.nf'
include { SAM_BAM_SORT_IDX } from '../../../process/samtools/convert-sort-index/main.nf'
include { QUAST } from '../../../process/quast/main.nf'

// include - workflow
include { SPADES_ABACAS } from '../../assembly/spades_abacas/main.nf'


workflow VIRAL_ASSEMBLY {
  // TODO: if usefull manage compression (.gz) of intermediate and final fasta ?

  take: 
  ch_trimmedReads   // channel: [ val(meta), [ path(reads, arity: 1..2) ] ] with at least 'meta.id', 'meta.host_kraken_db', 'meta.ref_genome' (optionally 'meta.spades_args')


  main:
  // host reads filtering
  if (!params.skip_host_filtering) {
    hostKrakenDb = ch_trimmedReads.map { it[0].host_kraken_db }
    KRAKEN2_HOST(ch_trimmedReads, hostKrakenDb)
    ch_cleanedReads = KRAKEN2_HOST.out.unclassified_reads_fastq
  } else {
    ch_cleanedReads = ch_trimmedReads
  }

  // assembly
  SPADES_ABACAS(ch_cleanedReads)
  scaffoldsAssembled = SPADES_ABACAS.out.scaffolds_res

  // mapping (with index building and bam sorting)
  if (!params.skip_bowtie2) {
    BOWTIE2_BUILD(scaffoldsAssembled)
    bwtIdx = BOWTIE2_BUILD.out.idx
    //// make reads and index channels in same order
    readsBowtieIndex = ch_cleanedReads.join(bwtIdx, remainder: true, by: 0)
    cleanedReads_chReordered = readsBowtieIndex.map { row -> tuple(row[0], row[1]) }
    bwtIdx_chReordered = readsBowtieIndex.map { row -> tuple(row[0], row[2]) }
    //TODO: when cleanedReads_chReordered.map { it[0].id }.collect() == bwtIdx_chReordered.map { it[0].id }.collect())
    BOWTIE2(cleanedReads_chReordered, bwtIdx_chReordered)
    rawSAM=BOWTIE2.out.sam
    SAM_BAM_SORT_IDX(rawSAM)
    sortedBAM = SAM_BAM_SORT_IDX.out.bam
  } else {
    sortedBAM = null
  }

  // assembly statistics (TODO: mark duplicates before)
  if (!params.skip_quast) {
    //// make inputRefGenome channels in same order as scaffoldsAssembled
    inputRefGenome = scaffoldsAssembled.map { it[0].ref_genome }
    QUAST(scaffoldsAssembled, sortedBAM, inputRefGenome)
    quastHtml = QUAST.out.html
  } else {
    quastHtml = null
  }

  emit:
  scaffoldsAssembled
  quastHtml
}


