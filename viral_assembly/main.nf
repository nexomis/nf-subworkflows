// include - process

include { FASTP_DEDUP } from '../../process/fastp/dedup/main.nf'
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

def nullifyEmptyFasta(fasta_file){
  def char_count = 0
  for (line in fasta_file.readLines()) {
    if (!line.startsWith(">")) {
      char_count += line.trim().length()
      if (char_count > 500) {
        return fasta_file
      }
    }
  }
  return null
}

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
      "do_abacas": ["String"],
      "keep_before_abacas": ["String"],
      "dedup": ["String"],
      "keep_before_dedup": ["String"]
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

  UNSPRING_READS(inputReads)
  UNSPRING_READS.out
  | map {
    new_meta = it[0].clone()
    new_meta.reads_id = "${new_meta.id}" // reads id definitoin
    new_meta.assembly_id = "${new_meta.reads_id}" // Assembly ID definition
    return [new_meta, it[1]]
  }
  | set { allFastq }

  // Start Remove Duplicates Before Assembly

  inputReadsUnclassified = ITERATIVE_UNCLASSIFIED_READS_EXTRACTION(allFastq, inputK2Index)

  inputReadsUnclassified.filter{ it[0].dedup == "yes" }
  | map { 
    new_meta = it[0].clone()
    new_meta.reads_id = "${new_meta.reads_id}_dedup"
    new_meta.assembly_id = "${new_meta.reads_id}"
    return [new_meta, it[1]]
  }
  | set { readsForDedup }

  dedupReads = FASTP_DEDUP(readsForDedup).reads

  inputReadsUnclassified.filter{ 
  (it[0].dedup == "no") || ((it[0].dedup == "yes") && (it[0].keep_before_dedup == "yes"))
  }
  | set { notDedupReadsToKeep }

  mergedForAssembly = notDedupReadsToKeep.concat(dedupReads)
  mergedForMapping = notDedupReadsToKeep.concat(readsForDedup) // dedup reads only for assembly

  // END Remove Duplicates Before Assembly

  mergedForAssembly 
  | flatMap { item ->
    def meta = item[0]
    def files = item[1]
    def assemblerList = meta.assembler
    assemblerList.collect { assemblerValue ->
        def new_meta = meta.clone()
        new_meta.assembler = assemblerValue
        [new_meta, files]
    }
  }
  | branch {
      spades: (it[0].assembler.startsWith("spades_"))
        def new_meta = it[0].clone()
        new_meta.args_spades = "${new_meta.assembler}".replace("spades_", "--")
        new_meta.assembly_id = "${new_meta.assembly_id}_${new_meta.assembler}"
        new_meta.label = "${new_meta.assembly_id}"
        return [new_meta, it[1]]
      no_assembly: true
    }
  | set { inputForAssembly }

  SPADES(inputForAssembly.spades)

  SPADES.out.scaffolds
  | map {
    def new_file = nullifyEmptyFasta(it[1])
    return [it[0].clone(), new_file] 
  }
  | filter { it[1] }
  | set { scaffolds }

  scaffolds
  | filter { it[0].do_abacas == "yes" }
  | map {[it[0].ref_id, it]}
  | filter { it[0] }
  | combine(faRefGenome.map {[it[0].id, it]}, by:0)
  | map {
      def new_meta = it[1][0].clone()
      new_meta.assembly_id = "${new_meta.assembly_id}_abacas"
      new_meta.label = "${new_meta.assembly_id}"
      return [it[0], [new_meta, it[1][1]], it[2]]
  }
  | set { joinInputForAbacas }

  // TODO SELECT BIGGER SCAFFOLD FOR NO_ABACAS ?? What about multi segment ?

  ABACAS(joinInputForAbacas.map {it[1]}, joinInputForAbacas.map {it[2]})

  ABACAS.out // process to extract the number of 
  | map {
    def new_meta = it[0]
    def new_file = nullifyEmptyFasta(it[1])
    return [new_meta, new_file]
  }
  | filter { it[1] } // remove empty abacas
  | concat(
      scaffolds.filter{ 
        it[0].keep_before_abacas == "yes" || it[0].do_abacas != "yes" 
      }
    )
  | set {finalScaffolds}

  // mapping (with index building and bam sorting)

  finalScaffolds
  | filter { it[0].realign == "yes" }
  | set { toRealignScaffolds }

  BOWTIE2_BUILD(toRealignScaffolds)
  BOWTIE2_BUILD.out.idx
  | set { bwtIdx }

  mergedForMapping.map {[it[0].id, it]}
  | combine(bwtIdx.map {[it[0].id, it]}, by:0) // merge with id (not dedup)
  | map {[[it[2][0], it[1][1]], it[2]]}
  | map {
      def new_meta = it[0][0]
      new_meta.label = "${new_meta.assembly_id}"
      return [[new_meta, it[0][1]], it[1]]
  }
  | set {joinInputforBt2}

  BOWTIE2(joinInputforBt2.map {it[0]},joinInputforBt2.map {it[1]})

  rawSAM = BOWTIE2.out
  SAM_BAM_SORT_IDX(rawSAM)
  sortedBAM = SAM_BAM_SORT_IDX.out.bam_bai

  emptyFile = EMPTY_FILE()

  refWithEmpty = faRefGenome.concat(emptyFile).map {[it[0].id, it]}
  bamWithEmpty = sortedBAM
  | map {["${it[0].assembly_id}", it]}
  | concat(emptyFile.map {[it[0].id, [it[0], it[1], it[1]]]})

  finalScaffolds
  | map { [(it[0].ref_id ? it[0].ref_id : "empty_file"), it] }
  | combine(refWithEmpty, by:0)
  | map {
    if (it[1][0].realign == "yes") {
      return ["${it[1][0].assembly_id}", it[1], it[2]]
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
      ids << item[0].assembly_id
      assemblies << item[1]
    }
    it[3].each{ item ->
      bams << item[1]
      bais << item[2]
    }
    meta = ["id": ids.join(","), label: it[0]]
    return [
      [meta, assemblies],
      it[2][0],
      [meta, bams, bais]
    ]
  }
  | set {inputForQuast}

  QUAST(inputForQuast.map{it[0]},inputForQuast.map{it[1]},inputForQuast.map{it[2]})

  emit:
  all_aln               = sortedBAM
  all_scaffolds         = finalScaffolds
  pre_abacas_scaffolds  = scaffolds
  post_abacas_scaffolds = ABACAS.out
  quast_dir             = QUAST.out.dir
  quast_html            = QUAST.out.html
  quast_tsv             = QUAST.out.tsv

}
