#!/usr/bin/env nextflow

// include - process
include { TRANSFERT_GFF } from '../../process/transfer-annot-viral/about-global-psa/main.nf'
include { BWA_INDEX } from '../../process/bwa/index/main.nf'
include { BOWTIE2_BUILD } from '../../process/bowtie2/bowtie2-build/main.nf'
include { BWA_MEM } from '../../process/bwa/mem/main.nf'
include { BOWTIE2 } from '../../process/bowtie2/mapping/main.nf'
include { SAM_BAM_SORT_IDX } from '../../process/samtools/convert-sort-index/main.nf'
include { PICARD_MARK_DUPLICATES } from '../../process/picard_tools/markduplicates/main.nf'
include { ABRA2 } from '../../process/abra2/main.nf'
include { SAV_CALL } from '../../process/sav_call/main.nf'
include { CALL_BATCH } from '../../process/sav_call/call_batch/main.nf'
include { checkMeta } from '../utils.nf'

workflow VIRAL_VARIANT {
  take:
  inputReads   // (meta, [fq_r1, {fq_r2}]): meta.id (sample_id), meta.batch_id, meta.rank_in_batch
  inputRef     // (meta_ref, [inRef_fa, {inRef_gff}]): meta.id (correspond to batch_id)
  inputAnnot   // (meta_annot, [inAnnot_fa, inAnnot_gff]): meta.id (correspond to batch_id)
  mapper       // available option: "bowtie2" or "bwa-mem"

  main:
  
  // Validate metadata structure
  def expectedMeta = [
    "batch_id": ["String"],
    "rank_in_batch": ["Integer"],
  ]
  inputReads.map { it -> checkMeta(it, expectedMeta) }
  inputRef.map { it -> checkMeta(it) }
  inputAnnot.map { it -> checkMeta(it) }

  // Validate input files
  inputReads.map { meta, file ->
    assert file.size() in 1..2 : "Invalid number of read files (${meta}, ${file})"
  }
  inputRef.map { meta, file ->
    assert file.size() in 1..2 : "Invalid number of reference files (${meta}, ${file})"
  }
  inputAnnot.map { meta, file ->
    assert file.size() == 2 : "Invalid number of annotation files (${meta}, ${file})"
  }

  // Validate batch concordance
  inputReads
  | map { it -> it[0].batch_id }
  | unique()
  | join(inputRef.map { it -> it[0].id }, remainder: true)
  | map { batch_id ->
      if (!batch_id) {
        error "Error: Different batches between reads and reference."
      }
    }

  // Validate annotation consistency
  inputRef
  | map { it -> [ it[0].id, it[1] ] }
  | join(inputAnnot.map { it -> [ it[0].id, it[1] ] }, remainder: true)
  | map { batch_id, ref, annot ->
    if (ref == null) {
      error "Error: reference undefined (${batch_id})"
    } else if (!( ((ref.size() == 1) && (annot != null)) ||
                  ((ref.size() == 2) && (annot == null)) )) {
        error "Error: annotation undefined or defined twice (at 'reference' and at 'annotation') (${batch_id})"
    }
  }

  // Transfer annotation from inputAnnot to inputRef when needed
  inputRef
  | map { it -> [ it[0].id, it ] }
  | branch { it -> 
    reannot: it[1][1].size() == 1
    with_own_annotation: it[1][1].size() == 2
  }
  | set { refAnnotStatus }

  // Process references needing annotation transfer
  refAnnotStatus.reannot
  | join(inputAnnot.map { it -> [ it[0].id, it ] })
  | map { _batch_id, ref, annot -> [ref, annot] }
  | set { toTransfer }

  TRANSFERT_GFF(toTransfer.map { it -> it[0] }, toTransfer.map { it -> it[1] })

  // Combine transferred annotations with references
  toTransfer
  | map { it -> [ it[0][0].id, it[0] ] }
  | join(TRANSFERT_GFF.out.transfered_gff.map { it -> [ it[0].id, it ] })
  | map { _batch_id, ref, gff -> [ ref[0], [ ref[1][0], gff[1] ] ] }
  | mix(refAnnotStatus.with_own_annotation.map { it -> it[1] })
  | set { refWithAnnot }

  // Perform read mapping based on selected mapper
  inputReads
  | map { it -> [ it[0].batch_id, it ] }
  | set { readsByBatch }

  if (mapper == "bowtie2") {
    BOWTIE2_BUILD(inputRef.map { it -> [ it[0], it[1][0] ] })
    
    readsByBatch
    | combine(BOWTIE2_BUILD.out.idx.map { it -> [ it[0].id, it ] }, by: 0)
    | map { _batch_id, reads, idx -> [reads, idx] }
    | set { toMap }

    BOWTIE2(toMap.map { it -> it[0] }, toMap.map { it -> it[1] })
    rawSam = BOWTIE2.out
  } else if (mapper == "bwa-mem") {
    BWA_INDEX(inputRef.map { it -> [ it[0], it[1][0] ] })
    
    readsByBatch
    | combine(BWA_INDEX.out.idx.map { it -> [ it[0].id, it ] }, by: 0)
    | map { _batch_id, reads, idx -> [reads, idx] }
    | set { toMap }

    BWA_MEM(toMap.map { it -> it[0] }, toMap.map { it -> it[1] })
    rawSam = BWA_MEM.out
  } else {
    error "Invalid mapper: '${mapper}' (supported: 'bowtie2', 'bwa-mem')"
  }

  // Convert SAM to BAM, sort and index
  SAM_BAM_SORT_IDX(rawSam)

  // Mark duplicates
  SAM_BAM_SORT_IDX.out.bam_bai
  | map { it -> [ it[0].batch_id, [ it[0], it[1] ] ] }
  | combine(inputRef.map { it -> [ it[0].id, [ it[0], it[1][0] ] ] }, by: 0)
  | map { _batch_id, bam, ref -> [bam, ref] }
  | set { toMarkDup }

  PICARD_MARK_DUPLICATES(toMarkDup.map { it -> it[0] }, toMarkDup.map { it -> it[1] })

  // Perform indel realignment
  PICARD_MARK_DUPLICATES.out.bam_bai
  | map { it -> [ it[0].batch_id, it ] }
  | combine(inputRef.map { it -> [ it[0].id, [ it[0], it[1][0] ] ] }, by: 0)
  | map { _batch_id, bam, ref -> [bam, ref] }
  | set { toRealign }

  ABRA2(toRealign.map { it -> it[0] }, toRealign.map { it -> it[1] })

  // Call variants per sample
  ABRA2.out.bam
  | map { it -> [ it[0].batch_id, [ it[0], it[1] ] ] }
  | combine(inputRef.map { it -> [ it[0].id, [ it[0], it[1][0] ] ] }, by: 0)
  | map { _batch_id, bam, ref -> [bam, ref] }
  | set { toCallVariants }

  SAV_CALL(toCallVariants.map { it -> it[0] }, toCallVariants.map { it -> it[1] })

  // Prepare batch variant calling inputs
  SAV_CALL.out.base_comp
  | map { meta, file -> [meta.batch_id, [meta.label ?: meta.id, file]] }
  | groupTuple()
  | join(
    SAV_CALL.out.snv_raw_indel
    | map { meta, file -> [meta.batch_id, [meta.label ?: meta.id, file]] }
    | groupTuple()
  )
  | map { batch_id, base_tuples, indel_tuples ->
    def meta = [
      id: batch_id,
      labels: base_tuples.collect { it -> it[0] }.join(',')
    ]
    [
      batch_id,
      [meta, base_tuples.collect { it -> it[1] }, indel_tuples.collect { it -> it[1] }]
    ]
  }
  | combine(inputRef.map { it -> [it[0].id, [it[0], it[1][0]]] }, by: 0)
  | combine(refWithAnnot.map { it -> [it[0].id, [it[0], it[1][1]]] }, by: 0)
  | map { _batch_id, batch_data, ref_data, gff_data ->
    [batch_data, ref_data, gff_data]
  }
  | set { toCallBatch }

  // Run batch variant calling
  CALL_BATCH(
    toCallBatch.map { it -> it[0] },
    toCallBatch.map { it -> it[1] },
    toCallBatch.map { it -> it[2] }
  )

  emit:
  // Per-sample outputs
  called_snv = SAV_CALL.out.snv_both
  called_snv_rev = SAV_CALL.out.snv_rev
  called_snv_fwd = SAV_CALL.out.snv_fwd
  called_indel = SAV_CALL.out.snv_raw_indel
  sav_base = SAV_CALL.out.base_comp
  
  // Batch-level outputs
  variants = CALL_BATCH.out.variants
  vcf = CALL_BATCH.out.vcf
  proteins = CALL_BATCH.out.proteins
  
  // Additional outputs
  transfered_gff = refWithAnnot.map { it -> it[1][1] }
  psa_align = TRANSFERT_GFF.out.psa
  flagstat = SAM_BAM_SORT_IDX.out.flagstat
  alignedBam = ABRA2.out.bam
}
