// process include
include { SPADES } from '../../../process/spades/main.nf'
include { ABACAS } from '../../../process/abacas/main.nf'


workflow SPADES_ABACAS {
  take: 
  ch_cleanedReads          // channel: [ val(meta), [ path(reads, arity: 1..2) ] ]

  main:
  // run spades
  SPADES(ch_cleanedReads)

  // run abacas
  if (!params.skip_abacas) {
    ch_spadesScaffolds = SPADES.out.scaffolds
    inputRefGenome = ch_spadesScaffolds.map { it[0].ref_genome }
    ABACAS(ch_spadesScaffolds, inputRefGenome)
    scaffolds = ABACAS.out.scaffolds
  } else {
    scaffolds = SPADES.out.scaffolds
  }
  emit:
  scaffolds_res = scaffolds
}
