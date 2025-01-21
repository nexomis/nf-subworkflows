#!/usr/bin/env nextflow
nextflow.preview.output = true


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
  inReads   // (meta, [fq_r1, {fq_r2}]): meta.id (sample_id), meta.batch_id, meta.rank_in_batch
  inRef     // (meta_ref, [inRef_fa, {inRef_gff}]): meta.id (correspond to batch_id)
  inAnnot   // (meta_annot, [inAnnot_fa, inAnnot_gff]): meta.id (correspond to batch_id)
  mapper    // available option: "bowtie2" or "bwa-mem"

  main:
  
  //// Validate input
  // validate metadata structure
  def expectedMeta = [
    "batch_id": ["String"],
    "rank_in_batch": ["Integer"],
  ]
  inReads.map { checkMeta(it, expectedMeta) }
  inRef.map { checkMeta(it) }
  inRef.map { checkMeta(it) } 
  // validate arity of paths of each inputs
  readsBatches = inReads.map { meta, file ->
    assert file.size() in 1..2 : "nb reads files not attempted (${meta}, ${file})"
  }
  refBatches = inRef.map { meta, file ->
    assert file.size() in 1..2 : "nb ref files not attempted (${meta}, ${file})"
  }
  annotBatches = inAnnot.map { meta, file ->
    assert file.size() == 2 : "nb annot files not attempted (${meta}, ${file})"
  }
  // validate concordance between batch of reads and batch of reference
  readsBatches = inReads.map {it[0].batch_id}
  refBatches = inRef.map { it[0].id }
  readsBatchesCount = readsBatches
                        .unique()
                        .count()
  refBatchesCount = refBatches.count()
  maxBatchesCount = readsBatchesCount
                      .join(refBatchesCount)
                      .max()
  unionBatchesCount = readsBatches
                        .unique()
                        .join(refBatches, remainder: true)
                        .count()
  unionBatchesCount
    .combine(maxBatchesCount)
    .map { union, max ->
      if (union != max) {
        error "Error: Different batches between reads and reference."
      }
    }

  // validate absence of duplicatation of ananotation (includde verification of concordance between batch of reference and potential batch of annotation)
  inRefIdTag = inRef.map { [ it[0].id, it[1] ] }
  inAnnotIdTag = inAnnot.map { [ it[0].id, it[1] ] }
  inRefIdTag
    .join(inAnnotIdTag, remainder: true)
    .map {batch_id, ref, annot ->
      if (ref == null) {
        error "Error: reference undefined (${batch_id})"
      } else if (!( ((ref.size() == 1) && (annot != null)) ||
                    ((ref.size() == 2) && (annot == null)) )) {
          error "Error: annotation undefined or defined twice (at 'reference' and at 'annotation') (${batch_id})"
      }
    }

  // Transfert annotation from inAnnot to inRef (when no inRef_gff) + genomic coord dict: g.inRefSNV -> g.inRefCoord
  // (TODO for last step of worflow: g.inRefCoord -> c.inRefCoord or coding -> coding to manage case where 1 locus is asigned to multiple coding element)
  refBatchTag = inRef.map { [ it[0].id, it ] }

  refBatchTagAnnotStatus = refBatchTag.branch {
                             reannot: it[1][1].size() == 1
                             with_own_annotation: it[1][1].size() == 2
                             other: error "Unexpected input encountered (error nb files given: '${it}')" }
  annotBatchTag = inAnnot.map { [ it[0].id, it ] }

  joinRefReannotAndAnnot = refBatchTagAnnotStatus.reannot.join(annotBatchTag, by: 0)

  TRANSFERT_GFF(joinRefReannotAndAnnot.map { it[1] },
                joinRefReannotAndAnnot.map { it[2] })

  transferedAnnot = TRANSFERT_GFF.out.transfered_gff
  refWithTransferedAnnot = joinRefReannotAndAnnot
                            .map { [ it[1][0].id, it[1] ] }
                            .join ( transferedAnnot.map { [ it[0].id, it ] }, by: 0)
                            .map{ [ it[1][0], [ it[1][1][0], it[2][1] ] ] }  // tuple (meta, [fa, gff])

  // mapping (with index building and bam sorting)
  readsBatchTag = inReads.map { [ it[0].batch_id, it ] }
  
  if ( mapper == "bowtie2" ) {
    BOWTIE2_BUILD(inRef.map { [ it[0], it[1][0] ] })
    bwtIdx = BOWTIE2_BUILD.out.idx

    bwtIdxBatchTag = bwtIdx.map { [ it[0].id, it ] }
    mergedInMapping = readsBatchTag.combine(bwtIdxBatchTag, by: 0)
    BOWTIE2(mergedInMapping.map { it[1] }, mergedInMapping.map { it[2] })
    rawSAM = BOWTIE2.out
  } else if ( mapper == "bwa-mem" ) {
    BWA_INDEX(inRef.map { [ it[0], it[1][0] ] })
    bwaIdx = BWA_INDEX.out.idx

    bwaIdxBatchTag = bwaIdx.map { [ it[0].id, it ] }
    mergedInMapping = readsBatchTag.combine(bwaIdxBatchTag, by: 0)
    BWA_MEM(mergedInMapping.map { it[1] }, mergedInMapping.map { it[2] })
    rawSAM = BWA_MEM.out
  } else {
    error "Invalid values for mapper: '${mapper}' (supported: 'bowtie2', 'bwa-mem')"
  }

  SAM_BAM_SORT_IDX(rawSAM)
  sortedBAM = SAM_BAM_SORT_IDX.out.bam_bai


  // remove/mark duplicates
  bamBatchTag = sortedBAM.map { [ it[0].batch_id, [ it[0], it[1] ] ] }  // excl. '.bai' (not useful for markduplcates)
  refBatchTag = inRef.map { [ it[0].id, [ it[0], it[1][0] ] ] }
  mergedInMarkDup = bamBatchTag.combine(refBatchTag, by: 0)

  PICARD_MARK_DUPLICATES( mergedInMarkDup.map{ it[1] }, mergedInMarkDup.map{ it[2] })
  deDupBAM = PICARD_MARK_DUPLICATES.out.bam_bai

  // indel realignment (include sort/index): Ideally, do it at batch level ?
  bamDeDupBatchTag = deDupBAM.map { [ it[0].batch_id, it ] }
  mergedInAbra2 = bamDeDupBatchTag.combine(refBatchTag, by: 0)

  ABRA2( mergedInAbra2.map{ it[1] }, mergedInAbra2.map{ it[2] })
  realignedBAM = ABRA2.out.bam


  // SAV_CALL per sample
  bamBatchTag = realignedBAM.map { [ it[0].batch_id, [ it[0], it[1] ] ] }  // excl. '.bai' (not useful for sav_call)
  refBatchTag = inRef.map { [ it[0].id, [ it[0], it[1][0] ] ] }
  mergedInSavCall = bamBatchTag.combine(refBatchTag, by: 0)
  
  SAV_CALL(mergedInSavCall.map{ it[1] }, mergedInSavCall.map{ it[2] })
  savCallSNV = SAV_CALL.out.snv_both
  savCallSNVrev = SAV_CALL.out.snv_rev
  savCallSNVfwd = SAV_CALL.out.snv_fwd
  savCallIndel = SAV_CALL.out.snv_raw_indel
  savCallbase = SAV_CALL.out.base_comp

  // Group SAV_CALL outputs by batch for CALL_BATCH
  baseFilesByBatch = savCallbase.map { meta, file -> 
    [meta.batch_id, [meta.label ?: meta.id, file]]
  }.groupTuple()

  indelFilesByBatch = savCallIndel.map { meta, file -> 
    [meta.batch_id, [meta.label ?: meta.id, file]]
  }.groupTuple()

  // Create batch metadata with sample labels
  batchMeta = baseFilesByBatch.join(indelFilesByBatch).map { batch_id, base_tuples, indel_tuples ->
    def labels = base_tuples.collect { it[0] }
    def base_files = base_tuples.collect { it[1] }
    def indel_files = indel_tuples.collect { it[1] }
    def meta = [id: batch_id, labels: labels.join(',')]
    [batch_id, meta, base_files, indel_files]
  }
  
  // Combine with reference and GFF
  gffBatchTag = refWithTransferedAnnot.map { meta, files -> 
    [meta.id, [meta, files[1]]]
  }

  refBatchTag = inRef.map { meta, files -> 
    [meta.id, [meta, files[0]]]
  }
  
  mergedInCallBatch = batchMeta
    .combine(refBatchTag, by: 0)
    .combine(gffBatchTag, by: 0)
    .map { batch_id, meta, base_files, indel_files, ref_tuple, gff_tuple ->
      tuple(
        tuple(meta, base_files, indel_files),
        tuple(ref_tuple[0], ref_tuple[1]),
        tuple(gff_tuple[0], gff_tuple[1])
      )
    }
  // Run batch variant calling
  CALL_BATCH(
    mergedInCallBatch.map { it[0] },
    mergedInCallBatch.map { it[1] },
    mergedInCallBatch.map { it[2] }
  )

  emit:   // conserve meta in emit ??
  sav_call_snv = savCallSNV
  sav_call_snv_rev = savCallSNVrev
  sav_call_snv_fwd = savCallSNVfwd
  sav_call_snv_indel = savCallIndel
  sav_call_snv_base = savCallbase
  variants = CALL_BATCH.out.variants
  vcf = CALL_BATCH.out.vcf
  proteins = CALL_BATCH.out.proteins
  transfered_gff = refWithTransferedAnnot.map { it[1][1] }
  psa_algn = TRANSFERT_GFF.out.psa
  flagstat = SAM_BAM_SORT_IDX.out.flagstat
  aln_bam = realignedBAM
}
