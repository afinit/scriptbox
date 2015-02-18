#!/usr/bin/env python
# functions for easy use in python


# returns a dictionary of the contigs found in the file passed
def get_contigs( infile ):
  contigs = {}

  # open file
  with open( infile, 'r' ) as fh:
    seq_id = line[0][1:]
    contigs[seq_id] = ''
    
    # loop through file
    for line in fh:
      line = line.split()
      # if line is a contig id, start new contig
      if line[0][0] == '>':
        seq_id = line[0][1:]
        contigs[seq_id] = ''
      # if line is not a contig id, append the sequence to the current id's sequence
      else:
        contigs[seq_id] += line[0]

  return contigs
