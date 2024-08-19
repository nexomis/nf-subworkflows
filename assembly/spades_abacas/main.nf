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

  // parse spades channels to abacas input : meta + scaffold_path
  spadesScaffolds = spadesOut[0].flatMap { tuple ->
    def (meta, output_dir) = tuple
    def spades_scaffolds_path = file("${output_dir}/${meta.id}.scaffolds.fa")
    return [[meta, spades_scaffolds_path]]
  }
  
  // run abacas
  ref_genome = ch_cleanedReads.map { it[0].ref_genome }
  abacasOut = ABACAS(spadesScaffolds, ref_genome, abacas_MUMmer_program, abacas_keep_on_output)

  //  parse abacas channels : meta + scaffold_path
  abacasScaffold = abacasOut[0].flatMap { tuple ->
    def (meta, output_dir) = tuple
    def abacas_scaffolds_path = file("${output_dir}/${meta.id}.fasta")
    return [[meta, abacas_scaffolds_path]]
  }

  emit:
  abacasScaffold
}


