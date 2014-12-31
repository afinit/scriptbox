#!/bin/python
# Purpose: Convert from amino acid sequence to each possible dna nucleotide sequence

import sys

# nucleotide to amino acid conversion
amino_dict = {
    'TGG' : 'W', 'TTT' : 'F', 'TTC' : 'F', 'TTA' : 'L', 'TTG' : 'L', 'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L',
    'ATT' : 'I', 'ATA' : 'I', 'ATC' : 'I', 'ATG' : 'M', 'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V',
    'TCT' : 'S', 'TCC' : 'S', 'TCA' : 'S', 'TCG' : 'S', 'AGT' : 'S', 'AGC' : 'S', 'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P',
    'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T', 'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A',
    'TAT' : 'Y', 'TAC' : 'Y', 'TAA' : '*', 'TAG' : '*', 'TGA' : '*', 'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q',
    'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG' : 'K', 'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E',
    'TGT' : 'C', 'TGC' : 'C', 'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R', 'AGA' : 'R', 'AGG' : 'R',
    'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G' }

def trans_amino( dna_seq ):
  # declare empty list
  amino_acid_list = []
  amino_acid_seq = ""
  amino_acid_seq_rev = ""

  for j in range( 0, len( dna_seq )/3 ):
    dna_trip = dna_seq[j*3:j*3+3]
    dna_trip_rev = dna_trip[::-1]
    amino_acid_seq += amino_dict[dna_trip]
    amino_acid_seq_rev += amino_dict[dna_trip_rev]
  
  amino_acid_list.append( amino_acid_seq )
  amino_acid_list.append( amino_acid_seq_rev )
  amino_acid_seq = ""
  amino_acid_seq_rev = ""

  for i in range( 1,3 ):
    for j in range( 0, len( dna_seq )/3-1 ):
      dna_trip = dna_seq[j*3+i:j*3+3+i]
      dna_trip_rev = dna_trip[::-1]
      amino_acid_seq += amino_dict[dna_trip]
      amino_acid_seq_rev += amino_dict[dna_trip_rev]
    
    amino_acid_list.insert( i, amino_acid_seq )
    amino_acid_list.insert( i+3, amino_acid_seq_rev )
    amino_acid_seq = ""
    amino_acid_seq_rev = ""

  return amino_acid_list

if len( sys.argv ) > 1 :
  print("|"+ sys.argv[1] +"|")
  amino_acid_list = trans_amino( sys.argv[1] )
  for i in range(6):
    print "Frame"+str(i)+": "+amino_acid_list[i]
