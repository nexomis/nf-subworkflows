include { MINIPROT } from '../../process/miniprot/main.nf'
include { checkMeta } from '../utils.nf'

process RENAME_PROT {
  container "python:3.10"

  label 'cpu_x1'
  label 'mem_2G'

  input:
  tuple val(meta), path(prot_fasta, arity: 1, stageAs: "inputs/*")

  output:
  tuple val(meta), path("${meta.id}.fa")

  script:
  """
  #!/usr/bin/env python3
  import re
  seq_dict = dict()
  seq_name_count = dict()
  with open("${prot_fasta}", "r") as in_file:
    for line in in_file.readlines():
      line = line.rstrip()
      if line.startswith(">"):
        match = re.match("${meta.regex_prot_name}", line)
        name = line.lstrip(">").split()[0]
        if match:
          if len(match.groups()) > 0:
            name = match.groups()[0]
        if name in seq_name_count.keys():
          seq_name_count[name] += 1
          if seq_name_count[name] == 2:
            seq_dict[name + "_1"] = seq_dict[name]
            del seq_dict[name]
          name = name + "_" + str(seq_name_count[name])
        else:
          seq_name_count[name] = 1
        seq_dict[name] = []
      else: 
        seq_dict[name].append(line)
  with open("${meta.id}.fa", "w") as ofile:
    for name, seqs in seq_dict.items():
      ofile.write(f">{name}\\n")
      for seq in seqs:
        ofile.write(seq+"\\n")
  """
}

process POST_PROCESS {
  container "quay.io/nexomis/pandas:2.1.4"

  label 'cpu_x1'
  label 'mem_2G'

  input:
  tuple val(meta_genome), path(genome_fasta, arity: 1, stageAs: "input/genome.fasta")
  tuple val(meta_annot), path(annot_file, arity: 1, stageAs: "input/annot.gff")

  output:
  tuple val(meta_genome), path("${meta_genome.label ?: meta_genome.id}.fasta"), emit: genome
  tuple val(meta_annot), path("${meta_annot.label ?: meta_annot.id}.gff"), emit: annot

  script:
  template "post_process.py"
}

workflow HANNOT {
    take:
    genomeFasta
    proteinFasta // check id proteinFasta ref_prot / id 

    main:

    def expectedMeta = [
      "proteome_id": ["String"],
      "revcomp": ["String"],
      "filter_annot": ["String"],
      "retain_only_annot": ["String"]
    ]

    def expectedMeta2 = [
      "regex_prot_name": ["String"]
    ]

    // Validate metadata
    genomeFasta.map { checkMeta(it, expectedMeta) }
    proteinFasta.map { checkMeta(it, expectedMeta2) }

    RENAME_PROT(proteinFasta)
    | set {renamedProt}

    genomeFasta.map{
      [it[0].proteome_id, it]
    }
    | combine(renamedProt.map{[it[0].id, it]}, by:0)
    | set {joinInputForMiniprot}

    MINIPROT(joinInputForMiniprot.map{it[1]}, joinInputForMiniprot.map{it[2]})
    | set {annotFile}

    genomeFasta.map{
      [it[0].id, it]
    }
    | combine(annotFile.map{[it[0].id, it]}, by:0)
    | set {joinPost}

    POST_PROCESS(joinPost.map{it[1]}, joinPost.map{it[2]})
    out_genome_file = POST_PROCESS.out.genome
    out_annot_file = POST_PROCESS.out.annot

    emit:
    genome = out_genome_file
    raw_annot = annotFile
    annot = out_annot_file

}
