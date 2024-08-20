// process include
include { KRAKEN2 } from '../../../process/kraken2/main.nf'

// workflow include
include { SPADES_ABACAS } from '../../assembly/spades_abacas/main.nf'

workflow VIRAL_ASSEMBLY {
  take: 
  ch_trimmedReads          // channel: [ val(meta), [ path(reads, arity: 1..2) ] ]
  kraken_db

  //now on : ext.conf abacas_MUMmer_program    // channel value: 'nucmed' or 'promer'
  //now using publish : abacas_keep_on_output    // channel value: 'scaffold', 'include_act_files', 'include_extended_act_files' or 'all'

  main:
  // run kraken2: host reads filtering
  if (!params.skip_host_filtering) {
    krakenOut = KRAKEN2(ch_trimmedReads, kraken_db)
    ch_cleanedReads = krakenOut.out.unclassified_reads_fastq
  } else {
    ch_cleanedReads = ch_trimmedReads
  }

  // run spades
  spadesOut = SPADES(ch_cleanedReads)

  // run abacas
  ch_spadesScaffolds = spadesOut.out.scaffolds
  ref_genome = ch_cleanedReads.map { it[0].ref_genome }
  abacasOut = ABACAS(spadesScaffolds, ref_genome, abacas_MUMmer_program, abacas_keep_on_output)

  emit:
  abacasOut.out.scaffold
}


