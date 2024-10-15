#!/usr/bin/env python3

# Options
FILTER_GFF = ${ meta_genome.filter_annot == "yes" ? "True" : "False"}                        # Option to filter the GFF
APPLY_REVCOMP = ${ meta_genome.revcomp == "yes" ? "True" : "False"}                          # Option to apply reverse complement
RETAIN_ANNOTATED_CONTIGS = ${ meta_genome.retain_only_annot == "yes" ? "True" : "False"}     # Option to retain only annotated contigs

# Input and Output Files
INPUT_GFF = "input/annot.gff"
INPUT_FASTA = "input/genome.fasta"
OUTPUT_GFF = "${meta_annot.label ?: meta_annot.id}.gff"
OUTPUT_FASTA = "${meta_genome.label ?: meta_genome.id}.fasta"

import pandas as pd
import re

# Complement dictionary
COMP = {"A":"T", "C":"G", "G":"C", "T":"A", "U":"A", "M":"K", "R":"Y",
      "W":"W", "S":"S", "Y":"R", "K":"M", "V":"B", "H":"D", "D":"H", "B":"V",
      "N":"N", "-":"-", ".":".", "!":"!", "?":"?"}

def read_gff(input_gff):
  gff3_col_names = ['seqid', 'source', 'type', 'start', 'end', 'score',
                    'strand', 'phase', 'attributes']
  df = pd.read_table(input_gff, comment="#", sep="\\t", names=gff3_col_names)
  df['ID'] = df['attributes'].str.extract(r'(?:ID=|Parent=)([^;]+)')
  return df

def filter_gff_annotations(df):
  mRNA_df = df[df['type'] == 'mRNA'].copy()
  # Extract attributes
  mRNA_df['Target'] = mRNA_df['attributes'].str.extract(r'Target=([^; ]+) \\d+ \\d+')
  mRNA_df['TargetStart'] = mRNA_df['attributes'].str.extract(r'Target=[^ ]+ (\\d+) \\d+').astype(int)
  mRNA_df['TargetEnd'] = mRNA_df['attributes'].str.extract(r'Target=[^ ]+ \\d+ (\\d+)').astype(int)
  mRNA_df['Identity'] = mRNA_df['attributes'].str.extract(r'Identity=([\\d.]+)').astype(float)
  # Calculate alignment length and score
  mRNA_df['AlignmentLength'] = (mRNA_df['TargetEnd'] - mRNA_df['TargetStart']).abs()
  mRNA_df['Score'] = mRNA_df['Identity'] * mRNA_df['AlignmentLength']
  # Resolve duplicates
  selected_ids_nested = mRNA_df.groupby('Target').apply(resolve_duplicates).reset_index(drop=True)
  # resolve_duplicates, how to get the flatten list without group info ?
  selected_ids = [item for sublist in selected_ids_nested.tolist() for item in sublist]
  filtered_df = df[df['ID'].isin(selected_ids)]
  return filtered_df

def resolve_duplicates(group):
  group = group.sort_values(by=['Score', 'Identity', 'AlignmentLength'], ascending=[False, False, False])
  selected = []
  while not group.empty:
      best_hit = group.iloc[0]
      selected.append(best_hit['ID'])
      # Remove overlapping hits
      group = group[(group['TargetEnd'] <= best_hit['TargetStart']) | (group['TargetStart'] >= best_hit['TargetEnd'])]
  return selected

def read_fasta(input_fasta):
  seq_dict = dict()
  with open(input_fasta, "r") as in_file:
      current_seq_name = None
      for line in in_file:
          line = line.rstrip()
          if line.startswith(">"):
              name = line.lstrip(">").split()[0]
              seq_dict[name] = ""
              current_seq_name = name
          else:
              seq_dict[current_seq_name] += line
  return seq_dict

def retain_annotated_contigs(seq_dict, df):
  annotated_contigs = df['seqid'].unique()
  seq_dict = {seqid: seq for seqid, seq in seq_dict.items() if seqid in annotated_contigs}
  return seq_dict

def determine_reverse_complement_contigs(df):
  is_minus = dict()
  for r in df.itertuples():
      if r.type == "CDS":
          if r.seqid in is_minus:
              is_minus[r.seqid] = is_minus[r.seqid] and (r.strand == "-")
          else:
              is_minus[r.seqid] = (r.strand == "-")
  return is_minus

def reverse_complement_sequences(seq_dict, is_minus):
  revcomp_len = dict()
  for name, do_revcomp in is_minus.items():
      if do_revcomp and name in seq_dict:
          seq = seq_dict[name]
          seq_rc = "".join([COMP.get(base.upper(), base.upper()) for base in seq[::-1]])
          seq_dict[name+"_revcomp"] = seq_rc
          del seq_dict[name]
          revcomp_len[name] = len(seq_rc)
  return [revcomp_len, seq_dict]

def adjust_positions(row, seq_lengths):
  if row['seqid'] in seq_lengths:
      seq_length = seq_lengths[row['seqid']]
      row['seqid'] = row['seqid']+"_revcomp"
      new_start = seq_length - row['end'] + 1
      new_end = seq_length - row['start'] + 1
      row['start'], row['end'] = new_start, new_end
      if row['strand'] == '+':
          row['strand'] = '-'
      elif row['strand'] == '-':
          row['strand'] = '+'
  return row

def write_fasta(seq_dict, output_fasta):
  with open(output_fasta, "w") as fh:
      for name, seq in seq_dict.items():
          fh.write(f">{name}\\n")
          for i in range(0, len(seq), 80):
              fh.write(seq[i:i+80].upper()+"\\n")

def write_gff(df, output_gff):
  with open(output_gff, "w") as fh:
    fh.write("##gff-version 3\\n")
  if not df.empty:
    df.to_csv(output_gff, sep='\\t', header=False, index=False, mode='a')

def main():
  # Read GFF file
  df = read_gff(INPUT_GFF)
  
  # Read FASTA file
  seq_dict = read_fasta(INPUT_FASTA)
  
  if not df.empty:
    # Filter GFF annotations if the option is enabled
    if FILTER_GFF:
        df = filter_gff_annotations(df)
  
    # Retain only annotated contigs if the option is enabled
    if RETAIN_ANNOTATED_CONTIGS:
        seq_dict = retain_annotated_contigs(seq_dict, df)

    # Apply reverse complement if the option is enabled
    if APPLY_REVCOMP:
        # Determine which contigs to reverse complement
        is_minus = determine_reverse_complement_contigs(df)
        # Reverse complement sequences and get their new lengths
        revcomp_len, seq_dict = reverse_complement_sequences(seq_dict, is_minus)
        # Adjust GFF positions
        df = df.apply(adjust_positions, axis=1, args=(revcomp_len,))
    else:
        revcomp_len = {}  # Empty dict if not reversing any sequences
  
  # Write the modified FASTA file
  write_fasta(seq_dict, OUTPUT_FASTA)
  
  # Write the adjusted GFF file
  write_gff(df, OUTPUT_GFF)

if __name__ == '__main__':
  main()