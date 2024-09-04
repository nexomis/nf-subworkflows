include { SPADES } from '../process/spades/main.nf'
include { ABACAS } from '../process/abacas/main.nf'

workflow {
  // Initialisation des canaux
  Channel
    .from([
      [ [ id: 'smpl_A', genome_ref: '/home/agoum/git/test/viral_assemnly/test_tmp/data/genome_ref.fasta', spades_args: '--rnaviral' ], 
        [ file('/home/agoum/git/test/viral_assemnly/test_tmp/data/smpl_A_R1.fastq'), file('/home/agoum/git/test/viral_assemnly/test_tmp/data/smpl_A_R2.fastq') ] ],
      [ [ id: 'smpl_B', genome_ref: '/home/agoum/git/test/viral_assemnly/test_tmp/data/genome_ref.fasta', spades_args: '--corona' ], 
        [ file('/home/agoum/git/test/viral_assemnly/test_tmp/data/smpl_B_R1.fastq') ] ]
    ])
    .set { ch_cleanedReads }

  Channel
    .from('nucmer')
    .set { abacas_MUMmer_program }
  
  Channel
    .from('scaffold')
    .set { abacas_keep_on_output }

  // Vérification du canal d'entrée
  // ch_cleanedReads.view { tuple ->
  //   println "ch_cleanedReads: ${tuple}"
  // }

  // Exécution du process SPADES
  scaffolds_raw = SPADES(ch_cleanedReads)

  // Vérification des sorties de SPADES
  // scaffolds_raw[0].view { tuple ->
  //  println "scaffolds_raw: ${tuple}"
  // }

  // Extraction du fichier de scaffolds
  ch_scaffolds = scaffolds_raw[0].flatMap { tuple ->
    def (meta, output_dir) = tuple
    def scaffolds_path = file("${output_dir}/${meta.id}.scaffolds.fa")
    return [[meta, scaffolds_path]]
  }

  // Vérification du canal ch_scaffolds
  //ch_scaffolds.view { tuple ->
  //  println "ch_scaffolds: ${tuple}"
  //}

  // Passage à ABACAS
  scaffold_abacas = ABACAS(ch_scaffolds, abacas_MUMmer_program, abacas_keep_on_output)

  emit:
  scaffold_abacas[0]
}
