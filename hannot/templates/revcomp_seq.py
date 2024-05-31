#!/usr/bin/env python3
import pandas as pd
COMP = {"A":"T", "C":"G", "G":"C", "T":"A", "U":"A", "M":"K", "R":"Y",
  "W":"W", "S":"S", "Y":"R", "K":"M", "V":"B", "H":"D", "D":"H", "B":"V",
  "N":"N", "-":"-", ".":".", "!":"!", "?":"?"}
gff3_col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']

seq_dict = dict()
with open("${genome_fasta}", "r") as in_file:
  for line in in_file.readlines():
    if line.startswith(">"):
      name = line.lstrip(">").split()[0]
      seq_dict[name] = ""
    else: 
      seq_dict[name] += line.rstrip()
is_only_minus = dict()

df = pd.read_table("${annot_file}", comment="#", sep="\\t",
  names=gff3_col_names)
is_minus = dict()
for r in df.itertuples():
  if r.type == "CDS":
    if r.seqid in is_minus.keys():
      is_minus[r.seqid] = is_minus[r.seqid] and r.strand == "-"
    else:
      is_minus[r.seqid] = r.strand == "-"

revcomp_len = dict()
for name, do_revcomp in is_minus.items():
  if do_revcomp:
    seq_dict[name] = "".join([COMP[x.upper()] for x in seq_dict[name][::-1]])
    revcomp_len[name] = len(seq_dict[name])

def adjust_positions(row, seq_lengths):
  if row['seqid'] in seq_lengths.keys():
    seq_length = seq_lengths[row['seqid']]
    new_start = seq_length - row['end'] + 1
    new_end = seq_length - row['start'] + 1
    row['start'], row['end'] = new_start, new_end
    if row['strand'] == '+':
      row['strand'] = '-'
    elif row['strand'] == '-':
      row['strand'] = '+'
  return row

df_adjusted = df.apply(adjust_positions, axis=1, args=(revcomp_len,))

with open("${name_genome}.fa", "w") as fh:
  for name, seq in seq_dict.items():
    fh.write(f">{name}\\n")
    for i in range(0, len(seq), 80):
      fh.write(seq[i:i+80].upper()+"\\n")

with open("${name_annot}.gff", "w") as fh:
  fh.write("##gff-version 3\\n")
df_adjusted.to_csv("${name_annot}.gff", sep='\\t', header=False, index=False,
  mode = "a")
