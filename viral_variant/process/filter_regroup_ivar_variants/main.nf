process FILTER_REGROUP_IVAR_VARIANTS {
  container "${params.biocontainers_registry ?: 'quay.io'}/biocontainers/pandas:2.2.1"   // 117 MiB !!

  input:
  tuple val(meta), path(ivar_tsv_files, arity: '1..*', stageAs: 'input/'), path(short_mpileup_cov_files, arity: '1..*', stageAs: 'input/')
  ////  BE CAREFUL: correspondance between 'tsv_file' and 'short_mpileup_file' is based on basename, not robust, prefer use list index (e.g: list ordered following sample[meta.rank_in_batch], or simply follong meta.id (but not ideal order in final file))


  output:
  tuple val(meta), path("*_corrected.tsv", type: 'file') , optional:false ,  emit: smpl_corrected_all_iSNVs
  tuple val(meta), path("*_filtered.tsv", type: 'file') , optional:false ,  emit: smpl_filtered_all_iSNVs
  tuple val(meta), path("${meta.label ?: meta.id}_summary_all_iSNVs.tsv", type: 'file') , optional:false ,  emit: batch_summary_all_iSNVs
  tuple val(meta), path("${meta.label ?: meta.id}_summary_all_iSNVs_light.tsv", type: 'file') , optional:false ,  emit: batch_summary_all_iSNVs_light
  tuple val(meta), path("${meta.label ?: meta.id}_summary_all_iSNVs_long_format.tsv", type: 'file') , optional:false ,  emit: batch_summary_all_iSNVs_long_frmt
  tuple val(meta), path("${meta.label ?: meta.id}_batchFiltered", type: 'dir') , optional:false ,  emit: smpl_batch_filtered_all_iSNVs

  script:
  min_dp = task.ext.min_dp ?: 30
  ref_dp_ratio_max = task.ext.ref_ratio_threshold ?: 0.9
  alt_dp_ratio_min = task.ext.alt_ratio_threshold ?: 0.002
  out_prefix = meta.label ?: meta.id

  template "filter_and_regroup_ivar_variants.py"

  stub:
  """
  #!/usr/bin/bash

  touch ${meta.label ?: meta.id}_summary_all_iSNVs.tsv stub_corrected.tsv stub_filtered.tsv
  mkdir ${meta.label ?: meta.id}_batchFiltered/

  """
}
