// process include
include { SPADES } from '../../../process/spades/main.nf'
include { ABACAS } from '../../../process/abacas/main.nf'


workflow SPADES_ABACAS {
  take: 
  ch_cleanedReads          // channel: [ val(meta), [ path(reads, arity: 1..2) ] ]
  abacas_MUMmer_program    // channel value: 'nucmed' or 'promer'
  abacas_keep_on_output    // channel value: 'scaffold', 'include_act_files', 'include_extended_act_files' or 'all'

  main:
  // run spades
  spadesOut = SPADES(ch_cleanedReads)

  // run abacas
  ch_spadesScaffolds = spadesOut.out.scaffolds
  ref_genome = ch_cleanedReads.map { it[0].ref_genome }
  abacasOut = ABACAS(ch_spadesScaffolds, ref_genome, abacas_MUMmer_program, abacas_keep_on_output)

  emit:
  abacasOut.out.scaffold
}


