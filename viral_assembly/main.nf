// include - process
include { BOWTIE2_BUILD } from '../../process/bowtie2/bowtie2-build/main.nf'
include { BOWTIE2 } from '../../process/bowtie2/mapping/main.nf'
include { SAM_BAM_SORT_IDX } from '../../process/samtools/convert-sort-index/main.nf'
include { QUAST } from '../../process/quast/main.nf'
include { GUNZIP } from '../../process/gunzip/main.nf'
include { SPADES } from '../../process/spades/main.nf'
include { ABACAS } from '../../process/abacas/main.nf'
include { EMPTY_FILE } from '../../process/empty_file/main.nf'
include { UNSPRING_READS } from '../unspring_reads/main.nf'
include { ITERATIVE_UNCLASSIFIED_READS_EXTRACTION } from '../iterative_unclassified_reads_extraction/main.nf'
include { checkMeta } from '../utils.nf'

workflow VIRAL_ASSEMBLY {

  take:
  inputReads
  inputK2Index
  inputRefGenome

  main:

  def expectedMeta = [
      "ref_id": ["String", "NullObject"],
      "class_db_ids": ["String[]", "ArrayList"],
      "class_tool": ["String"],
      "realign": ["String"],
      "do_abacas": ["String"]
  ]

  // Validate metadata
  inputReads.map { checkMeta(it, expectedMeta) }
  inputK2Index.map { checkMeta(it) }
  inputRefGenome.map { checkMeta(it) }

  sepRefGenome = inputRefGenome.branch {
    gzip: ["gz","gzip","z"].contains(it[1].getExtension().toLowerCase())
    other: true
  }

  GUNZIP(sepRefGenome.gzip)
  | concat(sepRefGenome.other)
  | set {faRefGenome}

  allFastq = UNSPRING_READS(inputReads)

  inputReadsUnclassified = ITERATIVE_UNCLASSIFIED_READS_EXTRACTION(allFastq, inputK2Index)

  inputReadsUnclassified 
  | flatMap { item ->
    def meta = item[0]
    def files = item[1]
    def assemblerList = meta.assembler
    assemblerList.collect { assemblerValue ->
        def newMeta = meta.clone()
        newMeta.assembler = assemblerValue
        [newMeta, files]
    }
  }
  | branch {
      spades: it[0].assembler.split(/:/)[0] == "spades"
        it[0].args_spades = "--" + it[0].assembler.split(/:/)[1]
        return it
      no_assembly: true
        return it
    }
  | set { inputForAssembly }

  SPADES(inputForAssembly.spades)

  SPADES.out.scaffolds
  | set { scaffolds }

  scaffolds
  | filter { it[0].do_abacas == "yes" }
  | map {[it[0].ref_id, it]}
  | filter { it[0] }
  | combine(faRefGenome.map {[it[0].id, it]}, by:0)
  | set { joinInputForAbacas }

  // TODO SELECT BIGGER SCAFFOLD FOR NO_ABACAS ?? What about multi segment ?

  ABACAS(joinInputForAbacas.map {it[1]}, joinInputForAbacas.map {it[2]})
  ABACAS.out
  | concat(scaffolds)
  | unique { it[0].id + it[0].assembler }
  | set {finalScaffolds}

  // mapping (with index building and bam sorting)

  finalScaffolds
  | filter { it[0].realign == "yes" }
  | set { toRealignScaffolds }

  BOWTIE2_BUILD(toRealignScaffolds)
  BOWTIE2_BUILD.out.idx
  | set { bwtIdx }

  inputReadsUnclassified.map {[it[0].id, it]}
  | combine(bwtIdx.map {[it[0].id, it]}, by:0)
  | map {[[it[2][0], it[1][1]], it[2]]}
  | set {joinInputforBt2}

  BOWTIE2(joinInputforBt2.map {it[0]}, joinInputforBt2.map {it[1]})

  rawSAM=BOWTIE2.out.sam
  SAM_BAM_SORT_IDX(rawSAM)
  SAM_BAM_SORT_IDX.out.bam
  | set { sortedBAM }

  emptyFile = EMPTY_FILE()
  //finalScaffolds.filter{it[0].realign != "yes"}

  refWithEmpty = faRefGenome.concat(emptyFile).map {[it[0].id, it]}
  bamWithEmpty = sortedBAM
  | map {["${it[0].id}:${it[0].assembler}", it]}
  | concat(emptyFile.map {[it[0].id, [it[0], it[1], it[1]]]})

  finalScaffolds
  | map { [(it[0].ref_id ? it[0].ref_id : "empty_file"), it] }
  | combine(refWithEmpty, by:0)
  | map {
    if (it[1][0].realign == "yes") {
      return ["${it[1][0].id}:${it[1][0].assembler}", it[1], it[2]]
    } else {
      return ["empty_file", it[1], it[2]]
    }
  }
  | combine(bamWithEmpty, by:0)
  | map {[it[1][0].id, it[1], it[2], it[3]]}
  | groupTuple()
  | map {
    def ids = []
    def assemblies = []
    def bams = []
    def bais = []
    it[1].each{ item ->
      ids << item[0].id + "_" + item[0].assembler
      assemblies << item[1]
    }
    it[3].each{ item ->
      bams << item[1]
      bais << item[2]
    }
    meta = ["id": ids.join(","), sample_id: it[0]]
    return [
      [meta, assemblies],
      it[2][0],
      [meta, bams, bais]
    ]
  }
  | set {inputForQuast}

  QUAST(inputForQuast.map{it[0]},inputForQuast.map{it[1]},inputForQuast.map{it[2]})

  emit:
  all_aln         = sortedBAM
  all_scaffolds   = finalScaffolds
  quast_dir       = QUAST.out.dir
  quast_html      = QUAST.out.html
  quast_tsv       = QUAST.out.tsv

}
