// Julien:
//  - Meta for each input of process (with at least meta.id!)
//  - Reduce the number of process calls by combining different processes/scripts! (The start and ending of tasks take more time!)
//  - Threshold on minimum proportion of majority position instead of threshold on minimum proportion of secondary abundance of nucleotides.
//  - Reference for SNP calling can be a sample!!!! So, in cases where ref_snp is a sample (=> different from ref_annotation), the definition of 'SNP' and 'coord' (always with the public reference) uses an independent process call instead of the same! Management of the same sample as reference for multiple batches?? Management of the same sample across multiple batches?
//  - Indel Warning: plot with ASCIIviz/SamView with a filter of concerned reads and zoom in on specific regions + score of enrichment in extremity read positions (requires samtools mpileup on filtered BAM?)
//  - Don't combine different nucleotide variations on the same codons.

// For now, variant calling is performed independently of strand (management of strand specificity complicates WIG parsing and output and presents too little biological interest: at worst, workflow can be run filtering input BAM by strand mapping, no?)
// What about uncovered regions in WIG???
// If indel realignement not inplémenté, evaluate the possibility of removing (as an option) the ends of reads in the bam (e.g. 10% in 5' and 3').

/////////////////////////
// TODO: Issue with the management of positions in the case of segmented viruses!! Add chromosome name at position at each step to manage this!! (PS: merging all segments with NNN isn't a proper solution...)
/////////////////////////

// For later:
//  - Indel realignment
//  - Management of strand specificity (with filtering of BAM at the first step)


//// VERSION 0:
// 1a. PAIRWISE_GLOBAL_ALIGN				      (0.fa_smpl,	0.fa_ref) 						      ->	coord + snp
// 1b. IGVTOOLS_COUNT						          (0.bam,		0.fa_ref) 						        ->	wig0
// 2.	 ADD_REF_POS_IN_WIG					        (1b.wig0,	1a.coord) 						        ->	wig1
// 3.	 ADD_STATS_IN_WIG					          (2.wig1)								              	->	wig2
// 4.	 IDENTIFY_HETEROGENEOUS_POS_IN_WIG 	(3.wig2)									              ->	heterPos
// 5.	 groovy_merge_pos_bySmpl_byBatch		(0.meta,	1a.snp,				4.heterPos)	    ->	varPosByBatch
// 6.  regroup_wigs_and_pos_by_batch  		(3.wig2,	5.varPosByBatch,	0.refFaGFF)	->	regroupWigRefPos
// 7.  SUMMARIZE_VARIANT_DETAILS_BY_BATCH (6.regroupWigRefPos)						        ->	wigSummaryBatch [+ VCF]
// 8.	 [SNF_EFF]
//// [VERSION 1:]
// 1a. PAIRWISE_GLOBAL_ALIGN				  (0.fa_smpl,	0.fa_ref) 						                ->	coord & snp
// 1b. IGVTOOLS_COUNT						      (0.bam,		0.fa_ref) 						                  ->	wig0
// 2.	 WIG_TO_SUMMARIZE_WIG_BY_BATCH  (0.meta,	1b.wig0,	1a.coord,	1a.snp,	0.refFaGFF)	-> wigSummaryBatch [+ VCF...]
// 3.  [SNF_EFF...]


//// new conventions:
// add stub process
// add meta for all input process
// refact input subworflow
// add validation of meta keys
// add worflow output & publish


// include - process
include { IGVTOOLS_COUNT } from '../../../process/igvtools/count/main.nf'


// PROCESS - PYTHON : nedd to be templated
process PAIRWISE_GLOBAL_ALIGN {
  input:
  tuple val(meta) path(sample_fa, ref_fa, arity: 2, stageAs: 'input/')  // sample_fa->alt_file   | ref_fa->ref_file

  output:
  tuple val(meta), path("${meta.id}_coords.txt", type: 'file') ,            optional:false ,  emit: coord
  tuple val(meta), path("${meta.id}_mutations.txt", type: 'file') ,         optional:false ,  emit: snp
  tuple val(meta), path("${meta.id}_globalAlgn.txt", type: 'file') ,        optional:true ,   emit: align
  tuple val(meta), path("${meta.id}_mutations_details.txt", type: 'file') , optional:true ,   emit: snp_details


  script:
  """
  python3 << OFE
  
  from Bio import SeqIO
  from Bio import Align

  ref_file=${ref_fa}
  alt_file=${sample_fa}
  smpl_name=${meta.id}
  verbose=True

  aligner = Align.PairwiseAligner()
  aligner.mode = "global"
  aligner.match_score = 5
  aligner.mismatch_score = -4
  aligner.open_gap_score = -10
  aligner.extend_gap_score = -0.5
  aligner.target_end_gap_score = -5
  aligner.query_end_gap_score = -5
  parsed_records = dict()
  parsed_records["ref"] = list(SeqIO.parse(ref_file, "fasta"))
  parsed_records["alt"] = list(SeqIO.parse(alt_file, "fasta"))
  records = dict()

  for key in parsed_records.keys():
    if len(parsed_records[key]) != 1:
      sys.exit("1 and only 1 contig per fasta file required")
    records[key] = parsed_records[key][0]

  if verbose:
    print("aligning %s with %s" % (records["ref"].id, records["alt"].id))

  alignments = aligner.align(records["ref"].seq, records["alt"].seq)
  best = alignments[0]
  if versboe:
    print(best)

  # Parse the alignment: coords + mutations
  alignment_array = best.__array__()
  ref = "".join(map(lambda x: x.decode("utf-8"), alignment_array[0]))
  alt = "".join(map(lambda x: x.decode("utf-8"), alignment_array[1]))

  coords = []
  mutations = []
  mutations_details = []
  if base0 :
    ref_pos = 0
    alt_pos = 0
  else :
    ref_pos = 1
    alt_pos = 1

  # TODO: add -1 if out of ref (5' or 3') using alt/ref_pos and length of ref/alt + add del/ins after end to prevent gaps at end ?
  track="start"
  for i in range(best.length):
    ref_char = ref[i]
    alt_char = alt[i]
    if ref_char != '-' and alt_char != '-':
      if ref_char != alt_char:
        mutations.append(ref_pos)
        mutations_details.append((ref_pos, alt_pos, 'Mismatch', ref_char, alt_char))
      coords.append((ref_pos, alt_pos))
      ref_pos += 1
      alt_pos += 1
      track="hit"
    elif ref_char == '-' and alt_char != '-':
      mutations.append(ref_pos)
      coords.append((ref_pos, alt_pos))
      if track == "start":
        mutations_details.append((ref_pos, alt_pos, 'Insertion_before_start', '-', alt_char))
      elif track == "ins":
        mutations_details.append((ref_pos, alt_pos, 'Insertion_extend', '-', alt_char))
        track="ins"
      else:
        mutations_details.append((ref_pos, alt_pos, 'Insertion_open', '-', alt_char))
        track="ins"
      alt_pos += 1
    elif ref_char != '-' and alt_char == '-':
      mutations.append(ref_pos)
      coords.append((ref_pos, alt_pos))
      if track == "start":
        mutations_details.append((ref_pos, alt_pos, 'Deletion_before_start', ref_char, '-'))
      elif track == "del":
        mutations_details.append((ref_pos, alt_pos, 'Deletion_extend', ref_char, '-'))
        track="del"
      else:
        mmutations_details.append((ref_pos, alt_pos, 'Deletion_open', ref_char, '-'))
        track="del"
      ref_pos += 1

  # Save results on files
  mutations_file = smpl_name + "_mutations.txt"
  with open(mutations_file, 'w') as f_mut:
      f_mut.write("Mutations positions:\n")
      for mutation in mutations:
          f_mut.write(f"{mutation}\n")

  coords_file = smpl_name + "_coords.txt"
  with open(coords_file, 'w') as f_coords:
      f_coords.write("Coordinates (ref_pos, alt_pos):\n")
      for coord in coords:
          f_coords.write(f"{c[0]},{c[1]}\n")

  if verbose:
    mutations_details_file = smpl_name + "_mutations_details.txt"
    with open(mutations_details_file, 'w') as f_details:
        f_details.write("Mutations details (ref_pos, alt_pos, type, ref_char, alt_char):\n")
        for mutation in mutations_details:
            f_details.write(f"{mutation[0]},{mutation[1]},{mutation[2]},{mutation[3]},{mutation[4]}\n")
    
    alignment_file = smpl_name + "_globalAlgn.txt"
    with open(alignment_file, 'w') as f_align:
      f_align.write("Best Alignment in plain text:\n")
      f_align.write(best.format('clustal'))
  
  EOF
  """
}

process ADD_REF_POS_IN_WIG {
  input:
  tuple val(meta), path(wigFile, coordFile, arity: 2, stageAs: 'input/')      // wig file, output of 'igvtools count' (unstranded !) AND coordoné between bam_ref (usually assembled genome) in 'PosInBamRef' column and reference genome (usually public reference) in 'PosInGeneralRef' column

  output:
  tuple val(meta), path("${meta.id}_w_ref.wig", type: 'file') , optional:false , emit: wig     // input wig enrichised by "Pos in Ref" column

  script:
  """
  python3 << OFE
  import csv

  with open('${coordFile}', 'r') as coord_f, open('${wigFile}', 'r') as wig_f:
    coords = {int(row['PosInBamRef']): int(row['PosInGlobalRef']) for row in csv.DictReader(coord_f)}
    output_lines = []

    for line in wig_f:
      if line.startswith('#') or line.startswith('track') or line.startswith('variableStep'):
        output_lines.append(line)
      else:
        fields = line.strip().split()
        pos_bam = int(fields[0])
        pos_ref = coords.get(pos_bam, 'NA')
        output_lines.append(f'{fields[0]}\t{pos_ref}\t' + '\t'.join(fields[1:]))
    
    with open('${meta.id}_w_ref.wig', 'w') as out_wig_f:
      for output_line in output_lines:
      out_wig_f.write(output_line + '\\n')
      
  EOF
  """
}

process ADD_STATS_IN_WIG {
  input:
  tuple val(meta), path(wigFile, arity: 1, stageAs: 'input/')      // wig file enrichised by "Pos in Ref" column

  output:
  tuple val(meta), path("${meta.id}_w_stats.wig", type: 'file') , optional:false , emit: wig     // input wig enrichised by coverage and proportions columns

  script:
  """
  python3 << OFE
  with open('${wigFile}', 'r') as wig_f:
    output_lines = []

    for line in wig_f:
      if line.startswith('#') or line.startswith('track') or line.startswith('variableStep'):
        output_lines.append(line)
      else:
        fields = line.strip().split()
        counts = [float(x) for x in fields[2:8]]
        cov = sum(counts)
        if cov == 0:
          stats = [0.0] * 6 + [cov]
        else:
          proportions = [count / cov for count in counts]
          stats = proportions + [cov]
        output_lines.append(line.strip() + '\\t' + '\\t'.join(map(str, stats)))

    with open('${meta.id}_w_stats.wig', 'w') as out_f:
      for output_line in output_lines:
        out_f.write(output_line + '\\n')
      
  EOF
  """
}

process IDENTIFY_HETEROGENEOUS_POS_IN_WIG {
    input:
    tuple val(meta), path(wigFile, arity: 1, stageAs: 'input/')      // wig file enrichised by "Pos in Ref" and statistics

    output:
    tuple val(meta), val(HETEROGENEOUS_POS) , optional:false , emit: list_pos  // list of heeterogeneous position ('Pos in Ref')

    script:
    """
    python3 << OFE
  # TODO: use column name instead of relative postition ?

  min_prop_variant = 0.1  # threshold for proportion of heterogenous status - ${task.ext.heterogn_prop_min ?: 0.1}
  min_cov = 30  # threshold for coverage of heterogenous status - ${task.ext.heterogn_cov_min ?: 30}

  with open('${wigFile}', 'r') as wig_f:
    positions = []
    for line in wig_f:
      if line.startswith('#') or line.startswith('track') or line.startswith('variableStep'):
        continue
      else:
        fields = line.strip().split()
        cov = float(fields[-1])
        if cov > min_cov :
          props = [float(x) for x in fields[-7:-1]]  # % A, C, G, T, INS, DEL
          second_abundance = sorted(props, reverse=True)[1]
          if second_abundance > min_prop_variant :
            pos_ref = fields[1]  
            positions.append(pos_ref)
      print(positions)    # return val instead of file this implies that this process process will be systematically executed, even with the '-resume'. TODO: in theory, it seems to me that the following steps of worflow should be compatible with option '-resume', but this need to be verified!
      
    EOF
    """
}

process SUMMARIZE_VARIANT_DETAILS_BY_BATCH {
  input:
  tuple val(meta), path(wig_files, arity: 1.., stageAs: 'input/')       // with meta.id=batch_id, and meta.sample_ids=lis of sample ids 
  path(ref_fa, ref_gff, arity: 2, stageAs: 'input/')
  val variantsPos           // no meta ???

  output:
  tuple val(meta.id), path("${meta.id}_variant_summary.tab", type: 'file') , optional:false , emit: tab

  script:
  """
  python3 << OFE

  from Bio import SeqIO
  from BCBio import GFF
  import pandas as pd

  ref_seq = SeqIO.read("${ref_fa}", "fasta")
  ref_gff = list(GFF.parse(open("${ref_gff}")))

  variant_positions = ${variantsPos}        ## for templating

  # summary df
  summary_df = pd.DataFrame(columns=["Ref_Pos", "Ref_Nucl", "Ref_ProtName", "Ref_ProtPos"])

  # ref inf. to summary df
  for pos in variant_positions:
    ref_nucl = ref_seq.seq[pos - 1]  # La séquence est indexée à partir de 0
    prot_name, prot_pos = "Intergenic", None
    for rec in ref_gff:
      for feature in rec.features:
        if feature.type == "gene":
          if feature.location.start <= pos <= feature.location.end:
          prot_name = feature.qualifiers.get('Name', ['Unknown'])[0]
          prot_pos = pos - feature.location.start + 1
    summary_df = summary_df.append({
      "Ref_Pos": pos,
      "Ref_Nucl": ref_nucl,
      "Ref_ProtName": prot_name,
      "Ref_ProtPos": prot_pos
    }, ignore_index=True)

  # parse wig: extract stat by smpl
  sample_stats = {}

  for sample_id, wig_file in zip(${meta.sample_ids}, ${wig_files}):
    # parse smpl wig
    bases_order = ['A', 'C', 'G', 'T', 'N', 'Del', 'Ins']
    with open(wig_file, 'r') as f:
      for line in f:
        if line.startswith('track'):
          continue
        elif line.startswith('variableStep'):
          continue
        #elif line.startswith('#Columns:'):
        #  cols = line.replace("^#Columns: ", "").split(',').strip()
        #  bases_order = [col.strip().replace("^Combined Strands ", "") for col in cols]
        else:
          cols = line.strip().split()
          pos = int(cols[0])
          if pos in variant_positions:
            bases_count = map(int, cols[2:7])  # A, C, G, T, N, DEL, INS
            consensus = max(zip(bases_order, bases_count), key=lambda x: x[1])[0]
            pA, pC, pG, pT, pN, pDel, pIns, cov = cols[8:14]
            pos_smpl = cols[1]
            
            sample_stats.setdefault(pos, {})[sample_id] = {
              'Cov': cov,
              'Consensus': consensus,
              'pA': pA, 'pC': pC, 'pG': pG, 'pT': pT,
              'pN': pN, 'pDel': pDel, 'pIns': pIns,
              'Pos_smpl': pos_smpl
            }

  # put smpl stats (from parsed wig) to summary_df
  for pos in variant_positions:
    if pos in sample_stats:
      for sample_id, stats in sample_stats[pos].items():
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_Cov'] = stats['Cov']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_Consensus'] = stats['Consensus']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_Pos'] = stats['pos_smpl']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_pA'] = stats['pA']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_pC'] = stats['pC']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_pG'] = stats['pG']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_pT'] = stats['pT']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_pN'] = stats['pN']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_pDel'] = stats['pDel']
        summary_df.loc[summary_df['Ref_Pos'] == pos, f'{sample_id}_pIns'] = stats['pIns']
    else :
      for sample_id in sample_ids:
        columns = [f'{sample_id}_Cov', f'{sample_id}_Consensus', f'{sample_id}_Pos'] + [f'{sample_id}_p{base}' for base in ['A', 'C', 'G', 'T', 'N', 'Del', 'Ins']]
        summary_df.loc[summary_df['Ref_Pos'] == pos, columns] = None
  # save summary_df
  summary_df.to_csv("${meta.id}_variant_summary.txt", sep="\t", index=False)

  EOF
  """
}


workflow VIRAL_VARIANT {

  take:
  inputs // (id, (meta, [reads_on_assembled_bam, {reads_on_assembled_bai}]), [assembled_fa, assembled_fai], [ref_genome_fa, ref_genome_gff]) with {xxx} = optional
  // TODO: meta.args_igvcount depend of params.variant_by_strand ?
  // TODO: meta.batch_id about samplesheet
  // !!!!!!!!! be carefull, by default fasta_sample, fasta_ref (gff_ref) and ref in bam are oriented in same strand !!!!!!!!!!!!

  main:
  
  // separate input tuple
  inputBam = inputs.map { [it[0], it[1]] }
  inputAssembledFa = inputs.filter { it[1] } | map { [it[0], it[2]] }
  inputRefGenomeFaGFF = inputs.filter { it[1] } | map { [it[0], it[3]] }

  // Validate metadata
  def expectedMeta = [
    "id": ["String", "NullObject"],
    "batch_id": ["String", "NullObject"]
  ]
  inputBam.map { checkMeta(it[1][0], expectedMeta) }

  // Position equivalence between assembled genome and general reference (coord) + SNP position (snp)
  inputAssembledFaOnly = inputAssembledFa.map { [it[0], it[1][0]] }
  inputRefGenomeFa = inputRefGenomeFaGFF.map { [it[0], it[1][0] }
  joinAssembledFaOnlyAndRefGenomeFa = inputAssembledFaOnly.join(inputRefGenomeFa, by:0)
  PAIRWISE_GLOBAL_ALIGN(joinAssembledFaOnlyAndRefGenomeFa)    // add meta in input channel !
  coord = PAIRWISE_GLOBAL_ALIGN.out.coord
  snp = PAIRWISE_GLOBAL_ALIGN.out.snp

  // Position abundances details (wig2) and heterogenous position (heterogenous_pos)
  joinInputBamAndAssembledFa = inputBam.join(inputAssembledFa, by:0)
  IGVTOOLS_COUNT(joinInputBamAndAssembledFa.map { it[1] }, joinInputBamAndAssembledFa.map { it[2] })
  wig0 = IGVTOOLS_COUNT.out.wig

  wig0withID = wig0.map { it[0].id, [it[0], it[1]] }
  joinWig0AndCoord = wig0withID.join(coord, by:0)
  ADD_REF_POS_IN_WIG(joinWig0AndCoord.map { it[1] }, joinWig0AndCoord.map { it[2] })
  wig1 = ADD_REF_POS_IN_WIG.out.wig
  
  ADD_STATS_IN_WIG(wig1)
  wig2 = ADD_STATS_IN_WIG.out.wig
  
  IDENTIFY_HETEROGENEOUS_POS_IN_WIG(wig2)
  heterogenousPosPath = IDENTIFY_HETEROGENEOUS_POS_IN_WIG.out.pos_list

  // Merge interest position by sample (snp + heterogenous_pos)
  heterogenousPos = heterogenousPosPath.map { meta, posFile ->
    def pos = posFile.text.readLines().collect { it as Integer }
    return [meta, pos] }
  heterogenousPosWithID = heterogenousPos.map { it[0].id, it[1] }

  joinSnpPosAndHeterogenousPos = snp.join(heterogenousPosWithID, by:0)
  variantsPosBySample = joinSnpPosAndHeterogenousPos.map { id, snpPos, heterogenousPos ->
    def mergedPos = (snpPos + heterogenousPos).unique().sort()
    return [id, mergedPos]}

  // Merge variants posisitons by batch
  id2Batch = inputs.map { [it[0], it[1][0].batch ] }
  joinVarPosBySmplAndId2Batch = variantsPosBySample.join(id2Batch, by: 0)
  varPosBySmplBatchTag = joinedChannel.map { it[2], it[1] }
  variantsPosByBatch = positionsByBatch.groupTuple(by: 0)
    .map { batch, listPos -> 
      def mergedPos = listPos.flatten().unique().sort()
      return tuple(batch, mergedPos) }
  
  // Extract variants details (variantsPosByBatch / wig2) and save it on file

  wig2ByBatch = wig2
    .map { meta, wig -> tuple(meta.batch_id, meta.id, wig) }
    .groupTuple(by: 0)
    .map { batch_id, list_of_meta_and_wigs ->
      def sample_ids = list_of_meta_and_wigs.collect { it[0] }
      def wig_paths = list_of_meta_and_wigs.collect { it[1] }
      def meta = [id: batch_id, sample_ids: sample_ids]
      tuple(val(meta), paths(wig_paths))    // path/paths when arity: 1 ?? need to be a list ?
    }
  //// if don't work:
  // wig2ByBatch = wig2
  //   .map { meta, wig -> tuple(meta.batch_id, meta.id, wig) }
  //   .groupTuple(by: 0)
  //   .map { batch_id, list_of_meta_and_wigs ->
  //     def sample_ids = list_of_meta_and_wigs.collect { it[0] }
  //     def wig_paths = list_of_meta_and_wigs.collect { it[1] }
  //     tuple(batch_id, tuple(sample_ids, wig_paths))
  //   }
  //   .map { batch_id, (sample_ids, wig_paths) ->
  //     def meta = [id: batch_id, sample_ids: sample_ids]
  //     tuple(val(meta), wig_paths)
  //   }

  sample2Batch = wig2
    .map { meta, wig -> tuple(meta.id, meta.batch_id) }
  
  batch2RefGenomeFaGFF = sample2Batch
    .join(inputRefGenomeFaGFF, by: 0)
    .map { it[0], it[2]}
    .unique()         // tuple ( val(batch_id), paths(ref_fa, ref_gff) )

  inputSummarizeVar = wig2ByBatch
    .map { meta, wig_paths -> val(meta.id), tuple(meta, wig_paths) }
    .join (variantsPosByBatch, by: 0)
    .join (batch2RefGenomeFaGFF, by: 0) // tuple( val(batch_id), tuple( val(meta), paths(wigs) ), path(ref_fa, ref_gff), val(variantsPosByBatch) ) 

  SUMMARIZE_VARIANT_DETAILS_BY_BATCH(inputSummarizeVar[1], inputSummarizeVar[2], inputSummarizeVar[3])

  // VCF about bam + PAIRWISE_GLOBAL_ALIGN.algn OR parse wig + PAIRWISE_GLOBAL_ALIGN.snp_details ??
  // ...


  emit:
  summary_nucl = SUMMARIZE_VARIANT_DETAILS_BY_BATCH.out.tab

}

