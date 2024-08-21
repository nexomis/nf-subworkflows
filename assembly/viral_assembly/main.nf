// include - process
include { KRAKEN2 } from '../../../process/kraken2/main.nf'
include { BOWTIE2_BUILD } from '../../../process/bowite2/bowtie2-build/main.nf'

// include - workflow
include { SPADES_ABACAS } from '../../assembly/spades_abacas/main.nf'


workflow VIRAL_ASSEMBLY {
  take: 
  ch_trimmedReads          // channel: [ val(meta), [ path(reads, arity: 1..2) ] ]
  kraken_db

  //now on : ext.conf abacas_MUMmer_program    // channel value: 'nucmed' or 'promer'
  //now using publish : abacas_keep_on_output    // channel value: 'scaffold', 'include_act_files', 'include_extended_act_files' or 'all'

  main:
  // host reads filtering
  if (!params.skip_host_filtering) {
    KRAKEN2(ch_trimmedReads, kraken_db)
    ch_cleanedReads = KRAKEN2.out.unclassified_reads_fastq
  } else {
    ch_cleanedReads = ch_trimmedReads
  }

  // assembly
  SPADES(ch_cleanedReads)

  // scaffolding
  ch_spadesScaffolds = SPADES.out.scaffolds
  inputRefGenome = ch_cleanedReads.map { it[0].ref_genome }
  ABACAS(spadesScaffolds, inputRefGenome)
  refGenome = ABACAS.out.scaffold

  // mapping (with build index ref)
  BOWTIE2_BUILD(refGenome)
  bwtIdx = BOWTIE2_BUILD.out.idx

  //// make channels of reads and index in same order
  readsBowtieIndex = ch_cleanedReads.join(bwtIdx, remainder: true, by: 0)
  cleanedReads_chReordered = readsBowtieIndex.map { row -> tuple(row[0], row[1]) }
  bwtIdx_chReordered = readsBowtieIndex.map { row -> tuple(row[0], row[2]) }
  //CHECK_CHANEL_ORDER(cleanedReads_chReordered.map { it[0].id }.collect(), bwtIdx_chReordered.map { it[0].id }.collect())

  BOWTIE2(cleanedReads_chReordered, bwtIdx_chReordered)
  rawSAM=BOWTIE2.out.raw_sam
  SAM_BAM_SORT_IDX(BOWTIE2.rawSAM)

  // assembly statistics (mark duplicates before !)
  sortedBAM = SAM_BAM_SORT_IDX.out.bam
  bai = SAM_BAM_SORT_IDX.out.bai  // same order than sortedBAM, so don't need to re-associated based on meta.id. Maybe is better to export bam and associated bai in same directories ?
  QUAST(refGenome, sortedBAM, bai?, ???)

  emit:
  refGenome
  QUAST.out.report
}


