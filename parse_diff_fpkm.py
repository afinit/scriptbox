#!/usr/bin/env python
# this will accept a gtf file and a gene_list file
# USAGE: parse_diff_fpkm.py -g gene_file.txt -d file.diff -s <sample> [-o <outfile>]

import os, sys
import argparse

def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='takes the genes listed in gene_file and prints the FPKM values of those genes matching sample to outfile from the diff file',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('-g,--gene_file', dest='gene_file', metavar='GENE_FILE', type=str, required=True, help='file containing the list of genes of interest')
  parser.add_argument('-d,--diff_file', dest='diff_file', metavar='DIFF_FILE', type=str, required=True, help='diff file from which to pull the FPKM values')
  parser.add_argument('-s,--sample', dest='sample', metavar='SAMPLE', type=str, required=True, help='sample from which to get the FPKM values')
  parser.add_argument('-o,--outfile', dest='outfile', metavar='OUTFILE', type=str, help='output file to which the FPKM values will be written, default value is the diff_file with sample appended and a .txt extension')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  gene_file = args.gene_file
  diff_file = args.diff_file
  sample = args.sample
  outfile = args.outfile

  # check value of outfile and handle if value wasn't specified
  if outfile is None:
    outfile = os.path.splitext( diff_file )[0] + '_' + sample + '.txt'

  # open gene_file and initialize gene_list
  gene_list_fh = open( gene_file, 'r' )
  gene_list = []

  # initialize line holding variable
  line = gene_list_fh.readline().rstrip()

  # loop through gene_file
  while line != '':
    gene_list.append( line )
    line = gene_list_fh.readline().rstrip()

  gene_list_fh.close()

  # open diff_file and get first line
  diff_fh = open( diff_file, 'r' )
  line = diff_fh.readline().split()

  # INITIALIZE
  # initialize sample_list, if 'sample' matches nothing this will show the user their options
  sample_list = set()
  no_match = True
  # initialize fpkm_list as dict with gene_list as keys with 0 as the default value
  fpkm_list = dict.fromkeys( gene_list, 0 )

  # loop through the diff_file and process data to find fpkm values of genes in gene_list.. if the genes are not found in the diff_file, the fpkm value will be 0
  while line != []:
    gene_name = line[2].rstrip('g')
    gene_name = gene_name.rstrip('.')

    if gene_name in gene_list:
      if line[4] == sample:
        no_match = False
        fpkm_list[gene_name] = line[7]
      elif line[5] == sample:
        no_match = False
        fpkm_list[gene_name] = line[8]
      else:
        sample_list.add(line[4])
        sample_list.add(line[5])

    line = diff_fh.readline().split()
   
  diff_fh.close()
  
  if no_match:
    print '{} not found in {}'.format( sample, diff_file )
    print 'Samples in {}:'.format( diff_file )
    for i in sorted(sample_list):
      print '\t{}'.format( i )
  
  else:
    with open( outfile, 'w' ) as out_fh:
      for i in fpkm_list:
        out_fh.write( i + '\t' + str(fpkm_list[i]) +'\n' )


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
