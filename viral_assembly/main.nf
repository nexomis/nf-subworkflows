// ============================================================================
// VIRAL ASSEMBLY WORKFLOW
// ============================================================================
// This workflow performs viral genome assembly from sequencing reads using:
// 1. Read preprocessing and deduplication
// 2. Anchor-based read selection
// 3. Taxonomic classification and filtering
// 4. De novo assembly with SPAdes
// 5. Post-assembly scaffolding (ABACAS) and annotation (HANNOT)
// 6. Quality assessment with QUAST

// Process imports - organized by functionality
include { FASTP_DEDUP } from '../../process/fastp/dedup/main.nf'
include { BOWTIE2_BUILD as BT2_BUILD; BOWTIE2_BUILD as BT2_BUILD_ANCHOR; BOWTIE2_BUILD as BT2_BUILD_REF } from '../../process/bowtie2/bowtie2-build/main.nf'
include { BOWTIE2 as BT2; BOWTIE2 as BT2_ANCHOR; BOWTIE2 as BT2_REF } from '../../process/bowtie2/mapping/main.nf'
include { SAM_BAM_SORT_IDX; SAM_BAM_SORT_IDX as SORT_ANCHOR; SAM_BAM_SORT_IDX as SORT_REF} from '../../process/samtools/convert-sort-index/main.nf'
include { EXTRACT_MAPPED_FASTQ } from '../../process/samtools/extract-mapped-fastq/main.nf'
include { CONCAT_FQ } from '../../process/concat_fq/main.nf'
include { QUAST } from '../../process/quast/main.nf'
include { GUNZIP as GUNZIP_ANCHOR; GUNZIP as GUNZIP_REF} from '../../process/gunzip/main.nf'
include { SPADES } from '../../process/spades/main.nf'
include { ABACAS } from '../../process/abacas/main.nf'
include { EMPTY_FILE } from '../../process/empty_file/main.nf'
include { UNSPRING_READS } from '../unspring_reads/main.nf'
include { HANNOT } from '../hannot/main.nf'
include { ITERATIVE_UNCLASSIFIED_READS_EXTRACTION } from '../iterative_unclassified_reads_extraction/main.nf'
include { checkMeta } from '../utils.nf'

def nullifyEmptyFasta(fasta_file){
  def char_count = 0
  def found_enough = false
  fasta_file.readLines().each { line ->
    if (!found_enough && !line.startsWith(">")) {
      char_count += line.trim().length()
      if (char_count > 500) {
        found_enough = true
      }
    }
  }
  return found_enough ? fasta_file : null
}

workflow VIRAL_ASSEMBLY {

  take:
  inputReads
  inputK2Index
  inputRefGenome
  anchorsFasta
  protFasta

  main:

  // ============================================================================
  // METADATA VALIDATION
  // ============================================================================
  // Define expected metadata fields and their types for input validation
  def expectedMeta = [
      "ref_id": ["String", "NullObject"],
      "anchor_id": ["String", "NullObject"],
      "class_db_ids": ["String[]", "ArrayList"],
      "class_tool": ["String"],
      "realign": ["String"],
      "do_abacas": ["String"],
      "keep_before_abacas": ["String"],
      "dedup": ["String"],
      "keep_before_dedup": ["String"]
  ]

  // Validate metadata for all input channels
  inputReads.map { it -> checkMeta(it, expectedMeta) }
  inputK2Index.map { it -> checkMeta(it) }
  inputRefGenome.map { it -> checkMeta(it) }
  anchorsFasta.map { it -> checkMeta(it) }
  protFasta.map { it -> checkMeta(it) }

  // ============================================================================
  // REFERENCE GENOME AND ANCHOR SEQUENCE PREPARATION
  // ============================================================================
  // Handle compressed reference genomes - decompress if needed
  sepRefGenome = inputRefGenome.branch { it -> 
    gzip: ["gz","gzip","z"].contains(it[1].getExtension().toLowerCase())
    other: true
  }

  GUNZIP_REF(sepRefGenome.gzip)
  | concat(sepRefGenome.other)
  | set {faRefGenome}

  // Handle compressed anchor sequences - decompress if needed
  sepAnchors = anchorsFasta.branch { it -> 
    gzip: ["gz","gzip","z"].contains(it[1].getExtension().toLowerCase())
    other: true
  }

  GUNZIP_ANCHOR(sepAnchors.gzip)
  | concat(sepAnchors.other)
  | set {faAnchors}

  // Build Bowtie2 index for anchor sequences (used for read selection)
  anchorsBt2Index = BT2_BUILD_ANCHOR(faAnchors)

  // Build Bowtie2 index for reference genomes (used for reference mapping)
  refBt2Index = BT2_BUILD_REF(faRefGenome)

  // ============================================================================
  // READ PROCESSING PRE-ANNOTATION
  // ============================================================================
  // Decompress reads and set up tracking IDs
  UNSPRING_READS(inputReads)
  UNSPRING_READS.out
  | map { it ->
    def new_meta = it[0].clone()
    new_meta.reads_id = "${new_meta.id}" // Unique identifier for this read set
    new_meta.assembly_id = "${new_meta.reads_id}" // Base assembly identifier
    return [new_meta, it[1]]
  }
  | set { allFastq }

  // ============================================================================
  // ANCHOR-BASED SELECTION
  // ============================================================================
  // Map reads to anchor sequences (if anchor_id is specified)
  // This helps select reads that are likely from the target viral genome
  allFastq
  | map { it ->[it[0].anchor_id, it]}
  | filter { it -> it[0]} // Only process samples with anchor_id
  | combine(anchorsBt2Index.map { it ->[it[0].id, it]}, by:0)
  | set {joinInputforBt2Anchors}

  BT2_ANCHOR(joinInputforBt2Anchors.map{ it -> it[1]}, joinInputforBt2Anchors.map{ it -> it[2]})
  | SORT_ANCHOR

  // Extract reads that mapped to anchors (target viral reads) and unmapped reads
  EXTRACT_MAPPED_FASTQ(SORT_ANCHOR.out.bam_bai.map{ it -> [it[0], it[1]]})

  // ============================================================================
  // TAXONOMIC CLASSIFICATION AND READ FILTERING
  // ============================================================================
  // Combine unmapped reads from anchor mapping with reads that had no anchor
  EXTRACT_MAPPED_FASTQ.out.unmapped
  | concat(allFastq.filter{ it -> !(it[0].anchor_id && it[0].anchor_id != "")})
  | set {inputForClassification}

  // Remove host/contaminant reads using iterative Kraken2 classification
  inputReadsUnclassified = ITERATIVE_UNCLASSIFIED_READS_EXTRACTION(inputForClassification, inputK2Index)
  
  // ============================================================================
  // MERGE ANCHOR-MAPPED READS WITH UNCLASSIFIED READS
  // ============================================================================
  // For samples with anchors: combine anchor-mapped reads with unclassified reads
  inputReadsUnclassified
  | filter { it -> it[0].anchor_id && it[0].anchor_id != ""}
  | map { it ->[it[0].id, it]}
  | combine(EXTRACT_MAPPED_FASTQ.out.mapped.map{ it -> [it[0].id, it]}, by:0)
  | set {joinInputForConcat}

  CONCAT_FQ(joinInputForConcat.map{ it -> it[1]}, joinInputForConcat.map{ it -> it[2]})
  | concat(inputReadsUnclassified.filter{ it -> !(it[0].anchor_id && it[0].anchor_id != "")})
  | set {inputReadsAfterSelection}

  // ============================================================================
  // OPTIONAL DEDUPLICATION
  // ============================================================================
  // Process reads for deduplication if requested
  inputReadsAfterSelection.filter{ it -> it[0].dedup == "yes" }
  | map { it -> 
    def new_meta = it[0].clone()
    new_meta.reads_id = "${new_meta.reads_id}_dedup"
    new_meta.assembly_id = "${new_meta.reads_id}"
    return [new_meta, it[1]]
  }
  | set { readsForDedup }

  dedupReads = FASTP_DEDUP(readsForDedup).reads

  // Keep reads that don't need deduplication or should be kept before dedup
  inputReadsAfterSelection.filter{ it -> 
  (it[0].dedup == "no") || ((it[0].dedup == "yes") && (it[0].keep_before_dedup == "yes"))
  }
  | set { notDedupReadsToKeep }

  inputReadsAfterSelection.filter{ it -> 
  (it[0].dedup == "no")
  }
  | set { notDedupReads }

  // Create separate channels for assembly (deduplicated) and mapping (original)
  mergedForAssembly = notDedupReadsToKeep.concat(dedupReads)
  mergedForMapping = notDedupReads.concat(readsForDedup) // Use original reads for mapping

  // ============================================================================
  // MAPPING TO THE REFERENCE
  // ============================================================================
  // Map reads to reference genome if ref_mapping is enabled and reference is available
  mergedForMapping
  | filter { it -> it[0].ref_id }
  | map { it -> [it[0].ref_id, it] }
  | combine(refBt2Index.map { it -> [it[0].id, it] }, by: 0)
  | map { it ->
      def ref_id = it[0]
      def reads_meta = it[1][0]
      def reads_files = it[1][1]
      def idx_meta = it[2][0]
      def idx_files = it[2][1]
      
      // Create new metadata for reference mapping
      def new_meta = reads_meta.clone()
      new_meta.label = "${new_meta.id}_ref_mapping"
      
      return [[new_meta, reads_files], [idx_meta, idx_files]]
  }
  | set { inputForRefMapping }

  // Map reads to reference genome
  BT2_REF(inputForRefMapping.map { it -> it[0] }, inputForRefMapping.map { it -> it[1] })

  // Convert SAM to sorted and indexed BAM
  SORT_REF(BT2_REF.out)
  refMappedBAM = SORT_REF.out.bam_bai

  mergedForMapping

  // ============================================================================
  // DE NOVO ASSEMBLY WITH SPADES
  // ============================================================================
  // Expand assembler parameters and prepare for multiple assembly strategies
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
  | branch { it -> 
      spades: (it[0].assembler.startsWith("spades_"))
        def new_meta = it[0].clone()
        new_meta.args_spades = "${new_meta.assembler}".replace("spades_", "--") // Convert to SPAdes arguments
        new_meta.assembly_id = "${new_meta.assembly_id}_${new_meta.assembler}"
        new_meta.label = "${new_meta.assembly_id}"
        return [new_meta, it[1]]
      no_assembly: true
    }
  | set { inputForAssembly }

  SPADES(inputForAssembly.spades)

  // Filter out empty or insufficient assemblies
  SPADES.out.scaffolds
  | map { it ->
    def new_file = nullifyEmptyFasta(it[1])
    return [it[0].clone(), new_file] 
  }
  | filter { it ->  it[1] } // Remove null (empty) assemblies
  | set { scaffolds }

  // ============================================================================
  // POST-ASSEMBLY SCAFFOLDING WITH ABACAS
  // ============================================================================
  // ABACAS orders and orients contigs against a reference genome
  scaffolds
  | filter { it ->  it[0].do_abacas == "yes" }
  | map { it ->[it[0].ref_id, it]}
  | filter { it ->  it[0] } // Only process if ref_id is specified
  | combine(faRefGenome.map { it ->[it[0].id, it]}, by:0)
  | map { it ->
      def new_meta = it[1][0].clone()
      new_meta.assembly_id = "${new_meta.assembly_id}_abacas"
      new_meta.label = "${new_meta.assembly_id}"
      return [it[0], [new_meta, it[1][1]], it[2]]
  }
  | set { joinInputForAbacas }

  // TODO: SELECT BIGGER SCAFFOLD FOR NO_ABACAS ?? What about multi segment ?

  ABACAS(joinInputForAbacas.map { it -> it[1]}, joinInputForAbacas.map { it ->  it[2]})

  // Filter out empty ABACAS results and combine with non-ABACAS scaffolds
  ABACAS.out
  | map { it ->
    def new_meta = it[0].clone()
    def new_file = nullifyEmptyFasta(it[1])
    return [new_meta, new_file]
  }
  | filter { it ->  it[1] } // Remove empty ABACAS results
  | concat(
      scaffolds.filter{ it -> 
        it[0].keep_before_abacas == "yes" || it[0].do_abacas != "yes" 
      }
    )
  | set {afterAbacasScaffolds}

  // ============================================================================
  // VIRAL ANNOTATION WITH HANNOT
  // ============================================================================
  // HANNOT performs homology-based annotation using protein databases
  HANNOT(
    afterAbacasScaffolds.filter{ it -> it[0].proteome_id != ""}.map{ it -> 
      def new_meta = it[0].clone()
      new_meta.assembly_id = "${new_meta.assembly_id}_hannot"
      new_meta.label = "${new_meta.assembly_id}"
      return [new_meta, it[1]]
    }, protFasta)

  // Filter out empty HANNOT results and combine with non-annotated scaffolds
  HANNOT.out.genome
  | map { it ->
    def new_meta = it[0].clone()
    def new_file = nullifyEmptyFasta(it[1])
    return [new_meta, new_file]
  }
  | filter { it ->  it[1] } // Remove empty annotation results
  | concat(
      afterAbacasScaffolds.filter{ it -> 
        it[0].keep_before_hannot == "yes" || it[0].proteome_id == "" 
      }
    )
  | set {finalScaffolds}

  // ============================================================================
  // READ MAPPING AND QUALITY ASSESSMENT
  // ============================================================================
  // Map reads back to final assemblies for quality assessment (if realign=yes)
  finalScaffolds
  | filter { it ->  it[0].realign == "yes" }
  | set { toRealignScaffolds }

  // Build Bowtie2 indices for final assemblies
  BT2_BUILD(toRealignScaffolds)
  BT2_BUILD.out.idx
  | set { bwtIdx }

  // Combine reads with corresponding assembly indices for mapping
  // Create mapping between original sample ID and reads
  mergedForMapping
  | map { it -> [it[0].id, it] }
  | set { readsById }

  // Create mapping between original sample ID and bowtie2 indices
  bwtIdx
  | map { it -> [it[0].id, it] }
  | set { indicesById }

  // Combine reads with their corresponding indices by original sample ID
  readsById
  | combine(indicesById, by: 0)
  | map { it ->
      def sample_id = it[0]
      def reads_meta = it[1][0]
      def reads_files = it[1][1]
      def idx_meta = it[2][0]
      def idx_files = it[2][1]
      
      // Create new metadata for mapping with assembly information
      def new_meta = reads_meta.clone()
      new_meta.assembly_id = idx_meta.assembly_id
      new_meta.label = idx_meta.assembly_id
      
      return [[new_meta, reads_files], [idx_meta, idx_files]]
  }
  | set { joinInputforBt2 }

  // Map reads to assemblies and convert to sorted BAM
  BT2(joinInputforBt2.map { it -> it[0] }, joinInputforBt2.map { it -> it[1] })

  rawSAM = BT2.out
  SAM_BAM_SORT_IDX(rawSAM)
  sortedBAM = SAM_BAM_SORT_IDX.out.bam_bai

  // ============================================================================
  // QUAST QUALITY ASSESSMENT PREPARATION
  // ============================================================================
  // Create empty file placeholder for samples without reference/alignment
  emptyFile = EMPTY_FILE()

  // Prepare reference genomes with empty file fallback
  refWithEmpty = faRefGenome.concat(emptyFile).map { it ->[it[0].id, it]}
  
  // Prepare BAM files with empty file fallback
  bamWithEmpty = sortedBAM
  | map { it ->["${it[0].assembly_id}", it]}
  | concat(emptyFile.map { it ->[it[0].id, [it[0], it[1], it[1]]]})

  // Prepare reference BAM files with empty file fallback
  refBamWithEmpty = refMappedBAM
  | map { it -> [it[0].id, [it[0], it[1], it[2]]] }
  | concat(emptyFile.map { it -> [it[0].id, [it[0], it[1], it[1]]] })

  // Complex channel operations to group assemblies, references, assembly BAMs, and reference BAMs for QUAST
  finalScaffolds
  | map { it -> [(it[0].ref_id ? it[0].ref_id : "empty_file"), it] }
  | combine(refWithEmpty, by:0)
  | map { it ->
    if (it[1][0].realign == "yes") {
      return ["${it[1][0].assembly_id}", it[1], it[2]]
    } else {
      return ["empty_file", it[1], it[2]]
    }
  }
  | combine(bamWithEmpty, by:0)
  | map { it ->
    def sample_id = it[1][0].id
    def ref_bam_key = it[1][0].ref_id ? sample_id : "empty_file"
    return [ref_bam_key, it[1], it[2], it[3]]
  }
  | combine(refBamWithEmpty, by:0)
  | map { it ->[it[1][0].id, it[1], it[2], it[3], it[4]]}
  | groupTuple()
  | map { it ->
    def ids = []
    def assemblies = []
    def assembly_bams = []
    def assembly_bais = []
    def ref_bams = []
    def ref_bais = []
    it[1].each{ item ->
      ids << item[0].assembly_id
      assemblies << item[1]
    }
    it[3].each{ item ->
      assembly_bams << item[1]
      assembly_bais << item[2]
    }
    def ref_bam = it[4][0][1]
    def ref_bai = it[4][0][2]
    def new_meta = ["id": ids.join(","), label: it[0]]
    return [
      [new_meta, assemblies],
      it[2][0],
      [new_meta, assembly_bams, assembly_bais],
      [new_meta, ref_bam, ref_bai]
    ]
  }
  | set {inputForQuast}

  // Run QUAST for assembly quality assessment
  QUAST(inputForQuast.map{ it -> it[0]},inputForQuast.map{ it -> it[1]},inputForQuast.map{ it -> it[2]},inputForQuast.map{ it -> it[3]})

  // ============================================================================
  // VIRAL ASSEMBLY OUTPUTS
  // ============================================================================
  // Declare outputs

  emit:
  all_aln               = sortedBAM
  ref_mapped_bam        = refMappedBAM
  cleaned_reads         = inputReadsAfterSelection
  anchored_reads        = EXTRACT_MAPPED_FASTQ.out.mapped
  unclassed_reads       = inputReadsUnclassified
  all_scaffolds         = finalScaffolds
  hannot_raw            = HANNOT.out.raw_annot
  hannot_filtered       = HANNOT.out.annot
  pre_abacas_scaffolds  = scaffolds
  post_abacas_scaffolds = ABACAS.out
  post_hannot_scaffolds = HANNOT.out.genome
  quast_dir             = QUAST.out.dir
  quast_html            = QUAST.out.html
  quast_tsv             = QUAST.out.tsv

}
