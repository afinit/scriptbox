#!/bin/python
# Purpose: Convert from amino acid sequence to each possible dna nucleotide sequence

import sys

# Amino Acid to dna conversion
dna_dict = {
  'W' : [ 'TGG' ],
  'F' : [ 'TTT', 'TTC' ],
  'L' : [ 'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG' ],
  'I' : [ 'ATT', 'ATC', 'ATA' ],
  'M' : [ 'ATG' ],
  'V' : [ 'GTT', 'GTC', 'GTA', 'GTG' ],
  'S' : [ 'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC' ],
  'P' : [ 'CCT', 'CCC', 'CCA', 'CCG' ],
  'T' : [ 'ACT', 'ACC', 'ACA', 'ACG' ],
  'A' : [ 'GCT', 'GCC', 'GCA', 'GCG' ],
  'Y' : [ 'TAT', 'TAC' ],
  '*' : [ 'TAA', 'TAG', 'TGA' ],
  'H' : [ 'CAT', 'CAC' ],
  'Q' : [ 'CAA', 'CAG' ],
  'N' : [ 'AAT', 'AAC' ],
  'K' : [ 'AAA', 'AAG' ],
  'D' : [ 'GAT', 'GAC' ],
  'E' : [ 'GAA', 'GAG' ],
  'C' : [ 'TGT', 'TGC' ],
  'R' : [ 'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG' ],
  'G' : [ 'GGT', 'GGC', 'GGA', 'GGG' ] }

def trans_dna( amino_acid_seq ):
  # initialize list with conversion of first amino acid
  dna = dna_dict[amino_acid_seq[0]]

  # Loop through the sequence of amino acids
  for i in range(1, len( amino_acid_seq ) ):
    amino_acid = dna_dict[amino_acid_seq[i]]
   
    # create copy of current list of amino_acids
    temp = dna[:]
    dna = []
    # loop through list of dna triplets corresponding to the current amino_acid
    for j in range(0, len( amino_acid )):
      dna += [s+amino_acid[j] for s in temp]

  return dna


if len( sys.argv ) > 1 :
  print("|"+ sys.argv[1] +"|")
  dna = trans_dna( sys.argv[1] )
  for s in dna:
    print s
