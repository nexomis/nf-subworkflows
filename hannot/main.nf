include { MINIPROT } from '../../process/miniprot/main.nf'

process RENAME_PROT {
  container "python:3.10"

  label 'cpu_x1'
  label 'mem_2G'

  input:
  tuple val(name), path(prot_fasta, arity: 1, stageAs: "inputs/*")

  output:
  tuple val(name), path("${name}.fa")

  script:
  """
  #!/usr/bin/env python3
  import re
  seq_dict = dict()
  seq_name_count = dict()
  with open("${prot_fasta}", "r") as in_file:
    for line in in_file.readlines():
      if line.startswith(">"):
        match = re.match("${params.regex_prot_name ?: '^>.*GN=([^ ]+).*\$'}", line)
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
  with open("${name}.fa", "w") as ofile:
    for name, seqs in seq_dict.items():
      ofile.write(f">{name}\\n")
      for seq in seqs:
        ofile.write(seq)
  """
}

process REVCOMP_SEQ {
  container "ghcr.io/nexomis/pandas:py3.11-2.0.3-1.0"

  label 'cpu_x1'
  label 'mem_2G'

  input:
  tuple val(name_genome), path(genome_fasta, arity: 1, stageAs: "input_genome/*")
  tuple val(name_annot), path(annot_file, arity: 1, stageAs: "input_annot/*")

  output:
  tuple val(name_genome), path("${name_genome}.fa"), emit: genome
  tuple val(name_annot), path("${name_annot}.gff"), emit: annot

  script:
  template "revcomp_seq.py"
}

workflow HANNOT {
    take:
    genomeFasta
    proteinFasta

    main:

    // take basename if not given

    namedProteinFasta = proteinFasta.map {
      if (it instanceof Path) {
        return [it.getBaseName(), it]
      } else {
        return it
      }
    }

    namedGenomeFasta = genomeFasta.map {
      if (it instanceof Path) {
        return [it.getBaseName(), it]
      } else {
        return it
      }
    }

    RENAME_PROT(namedProteinFasta)
    | set {renamedProt}

    MINIPROT(namedGenomeFasta, namedProteinFasta)
    | set {annotFile}

    if (params.revcomp_from_annot) {
      REVCOMP_SEQ(namedGenomeFasta, annotFile)
      out_genome_file = REVCOMP_SEQ.out.genome
      out_annot_file = REVCOMP_SEQ.out.annot
    } else {
      out_genome_file = namedGenomeFasta
      out_annot_file = annotFile
    }

    emit:
    genome = out_genome_file
    annot = out_annot_file

}
