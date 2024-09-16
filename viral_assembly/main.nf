// include - process
include { KRAKEN2 as KRAKEN2_HOST1; KRAKEN2 as KRAKEN2_HOST2; KRAKEN2 as KRAKEN2_HOST3 } from '../../process/kraken2/main.nf'
include { GZ as GZ1; GZ as GZ2; GZ as GZ3 } from '../../process/gz/main.nf'
include { BOWTIE2_BUILD } from '../../process/bowtie2/bowtie2-build/main.nf'
include { BOWTIE2 } from '../../process/bowtie2/mapping/main.nf'
include { SAM_BAM_SORT_IDX } from '../../process/samtools/convert-sort-index/main.nf'
include { QUAST } from '../../process/quast/main.nf'
include { SPADES } from '../../process/spades/main.nf'
include { ABACAS } from '../../process/abacas/main.nf'
include { EMPTY_FILE } from '../../process/empty_file/main.nf'
include { SPRING_DECOMPRESS } from '../../process/spring/decompress/main.nf'
include { checkMeta } from '../utils.nf'

workflow VIRAL_ASSEMBLY {

  take:
  inputReads
  inputK2Index
  inputRefGenome

  main:

  def expectedMeta = [
      "ref_id": ["String", "NullObject"],
      "k2_ids": ["String[]", "ArrayList"]
  ]

  // Validate metadata
  inputReads.map { checkMeta(it, expectedMeta) }
  inputK2Index.map { checkMeta(it) }
  inputRefGenome.map { checkMeta(it) }

  inputReads
  | branch {
    spring: it[0].read_type == "spring"
    fastq: true
  }
  | set {reads}

  reads.spring  
  | SPRING_DECOMPRESS
  | set {unspringReads}

  reads.fastq
  | concat(unspringReads)
  | set {allFastq}

  allFastq.map { [it[0].k2_ids.size() >= 1 ? it[0].k2_ids[0] : "", it] }
  | filter { it[0] != "" }
  | combine(inputK2Index.map{[it[0].id, it]}, by:0)
  | set {joinInputForK2i1}

  KRAKEN2_HOST1(joinInputForK2i1.map { it[1] }, joinInputForK2i1.map { it[2] })
  KRAKEN2_HOST1.out.unclassified_reads_fastq
  | GZ1
  | concat(inputReads)
  | unique { it[0].id }
  | set {inputReadsFromK1}

  inputReadsFromK1.map { [it[0].k2_ids.size() >= 2 ? it[0].k2_ids[1] : "", it] }
  | filter { it[0] != "" }
  | combine(inputK2Index.map{[it[0].id, it]}, by:0)
  | set {joinInputForK2i2}

  KRAKEN2_HOST2(joinInputForK2i2.map { it[1] }, joinInputForK2i2.map { it[2] })
  KRAKEN2_HOST2.out.unclassified_reads_fastq
  | GZ2
  | concat(inputReadsFromK1)
  | unique { it[0].id }
  | set {inputReadsFromK2}

  inputReadsFromK2.map { [it[0].k2_ids.size() >= 3 ? it[0].k2_ids[2] : "", it] }
  | filter { it[0] != "" }
  | combine(inputK2Index.map{[it[0].id, it]}, by:0)
  | set {joinInputForK2i3}

  KRAKEN2_HOST3(joinInputForK2i3.map { it[1] }, joinInputForK2i3.map { it[2] })
  KRAKEN2_HOST3.out.unclassified_reads_fastq
  | GZ3
  | concat(inputReadsFromK2)
  | unique { it[0].id }
  | set {inputReadsFromK3}

  inputReadsFromK3 
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
  | map {[it[0].ref_id, it]}
  | filter { it[0] }
  | combine(inputRefGenome.map {[it[0].id, it]}, by:0)
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

  inputReadsFromK3.map {[it[0].id, it]}
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

  refWithEmpty = inputRefGenome.concat(emptyFile).map {[it[0].id, it]}
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
