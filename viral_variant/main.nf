#!/usr/bin/env nextflow
nextflow.preview.output = true


// include - process
include { TRANSFERT_GFF } from '../../process/transfer-annot-viral/about-global-psa/main.nf'
include { BOWTIE2_BUILD } from '../../process/bowtie2/bowtie2-build/main.nf'
include { BOWTIE2 } from '../../process/bowtie2/mapping/main.nf'
include { SAM_BAM_SORT_IDX } from '../../process/samtools/convert-sort-index/main.nf'
include { ABRA2 } from '../../process/abra2/main.nf'
include { IVAR_VARIANTS_ALL } from '../../process/ivar/variants/main.nf'
include { FILTER_REGROUP_IVAR_VARIANTS } from './process/filter_regroup_ivar_variants/main.nf'
include { checkMeta } from '../utils.nf'




workflow VIRAL_VARIANT {
  take:
  inReads   // (meta, [fq_r1, {fq_r2}]): meta.id (sample_id), meta.batch_id, meta.rank_in_batch
  inRef     // (meta_ref, [inRef_fa, {inRef_gff}]): meta.id (correspond to batch_id)
  inAnnot  // (meta_annot, [inAnnot_fa, inAnnot_gff]): meta.id (correspond to batch_id)


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

  // for upcoming version: reads correction (already performed in spades?): really useful ?

  // mapping (with index building and bam sorting)
  BOWTIE2_BUILD(inRef.map { [ it[0], it[1][0] ] })
  bwtIdx = BOWTIE2_BUILD.out.idx

  readsBatchTag = inReads.map { [ it[0].batch_id, it ] }
  bwtIdxBatchTag = bwtIdx.map { [ it[0].id, it ] }
  mergedInMapping = readsBatchTag.combine(bwtIdxBatchTag, by: 0)
  BOWTIE2(mergedInMapping.map { it[1] }, mergedInMapping.map { it[2] })

  rawSAM = BOWTIE2.out
  SAM_BAM_SORT_IDX(rawSAM)
  sortedBAM = SAM_BAM_SORT_IDX.out.bam_bai


  // for upcoming version: mark duplicate


  // indel realignment (include sort/index): Ideally, do it at batch level ?
  bamBatchTag = sortedBAM.map { [ it[0].batch_id, it ] }
  refBatchTag = inRef.map { [ it[0].id, [ it[0], it[1][0] ] ] }
  mergedInAbra2 = bamBatchTag.combine(refBatchTag, by: 0)

  ABRA2( mergedInAbra2.map{ it[1] }, mergedInAbra2.map{ it[2] })
  realignedBAM = ABRA2.out.bam

  // for upcoming version: ivar consensus - from specific mpileup to set specific parameters (quality, depth, ...)


  // unfiltered mpileup + unfiltered ivar
  completeRefBatchTag = refBatchTagAnnotStatus.with_own_annotation
                          .mix( refWithTransferedAnnot
                            .map { [ it[0].id, it ] } )

  bamWithRefBatchTag = realignedBAM
                         .map{ [ it[0].batch_id, [ it[0], [it[1], it[2]] ] ] }
                         .combine(completeRefBatchTag, by: 0)

  IVAR_VARIANTS_ALL( bamWithRefBatchTag.map { it[1] },
                         bamWithRefBatchTag.map { it[2] } )

  ivarRawBySmpl = IVAR_VARIANTS_ALL.out.raw_iSNV_tsv
  mpileupCovBySmpl = IVAR_VARIANTS_ALL.out.mpileup_cov

  // filter ivar_complete to keep only interest pos and write summary file by batch (following order designed by 'rank_in_batch' key)

  ivarRawFilesByBatch = ivarRawBySmpl
                          .map { [it[0].rank_in_batch as Integer, it[0].batch_id, it[1] ] }
                          .toSortedList{ a, b -> a[0] <=> b[0] }
                          .flatMap()
                          .map { [ it[1], it[2] ] }
                          .groupTuple()

  mpileupCovFilesByBatch = mpileupCovBySmpl
                             .map { [it[0].rank_in_batch as Integer, it[0].batch_id, it[1] ] }
                             .toSortedList{ a, b -> a[0] <=> b[0] }
                             .flatMap()
                             .map { [ it[1], it[2] ] }
                             .groupTuple()

  joinIvarAndMpileupByBatch = ivarRawFilesByBatch
                                .join(mpileupCovFilesByBatch)
                                .map{ [ [id:it[0]], it[1], it[2] ] }

  FILTER_REGROUP_IVAR_VARIANTS( joinIvarAndMpileupByBatch )

  // add 'IVAR_VARIANTS_ALL.out.mpileup_cov' as input of 'FILTER_REGROUP_IVAR_VARIANTS' to fill the blanks REF_DP at unvaribale position/sample
  summaryVarByBatch = FILTER_REGROUP_IVAR_VARIANTS.out.batch_summary_all_iSNVs
  summaryVarByBatchLight = FILTER_REGROUP_IVAR_VARIANTS.out.batch_summary_all_iSNVs_light
  summaryVarByBatchLgFrmt = FILTER_REGROUP_IVAR_VARIANTS.out.batch_summary_all_iSNVs_long_frmt


  // for upcoming version: in option (at batch level when 'inAnnot' !) edit coord in batch_summary_all_iSNVs (and individual file ?): g.inAnnot -> g.inRefCoord and g.inRefCoord -> c.inRefCoord
  // (with management of locus assigned to multiple genes!!! better genes_ref -> genes_alt???)


  // for upcoming version: plots QC and results (cf. entourage: 'https://codeberg.org/CENMIG/entourage/src/branch/v1.0/example%20results')
  // (% and nb reads mapped, alignement score, % bases assigned, cov of each smpl by batch, % base covered, % variable base uncovered on at least one smpl (by batch), indel warning)


  emit:   // conserve meta in emit ??
  summary_var_by_batch = summaryVarByBatch
  summary_var_by_batch_light = summaryVarByBatchLight
  summary_var_by_batch_long_frmt = summaryVarByBatchLgFrmt
  var_batch_filtered = FILTER_REGROUP_IVAR_VARIANTS.out.smpl_batch_filtered_all_iSNVs             //
  var_by_smpl_filtered = FILTER_REGROUP_IVAR_VARIANTS.out.smpl_filtered_all_iSNVs
  var_by_smpl_corrected = FILTER_REGROUP_IVAR_VARIANTS.out.smpl_corrected_all_iSNVs        
  subset_mpileup_cov_by_smpl = mpileupCovBySmpl     //
  transfered_gff = refWithTransferedAnnot.map { it[1][1] }
  psa_algn = TRANSFERT_GFF.out.psa                         //
  psa_genomic_coords = TRANSFERT_GFF.out.genomic_coords    //
  flagstat = SAM_BAM_SORT_IDX.out.flagstat
  aln_bam = realignedBAM
}

