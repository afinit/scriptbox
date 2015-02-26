#!/usr/bin/env python
# this will accept a gtf file and a gene_list file
# USAGE: parse_gtf.py gene_file.txt transcripts.gtf

import os, sys
import argparse

def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='takes the genes listed in gene_file and prints the FPKM values of those genes to outfile from the gtf file',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('-g,--gene_file', dest='gene_file', metavar='GENE_FILE', type=str, required=True, help='file containing the list of genes of interest')
  parser.add_argument('-f,--gtf_file', dest='gtf_file', metavar='GTF_FILE', type=str, required=True, help='gtf file from which to pull the FPKM values')
  parser.add_argument('-o,--outfile', dest='outfile', metavar='OUTFILE', type=str, help='output file to which the FPKM values will be written, default value is the gtf_file with a .txt extension')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  gene_file = args.gene_file
  gtf_file = args.gtf_file
  outfile = args.outfile

  # check value of outfile and handle if value wasn't specified
  if outfile is None:
    outfile = os.path.splitext( gtf_file )[0] + '.txt'

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

  # initialize fpkm_list as dict with gene_list as keys with 0 as the default value
  fpkm_list = dict.fromkeys( gene_list, 0 )

  # open getf_file and get first line
  gtf_fh = open( gtf_file, 'r' )
  line = gtf_fh.readline().split()

  # loop through the gtf_file and process data to find fpkm values of genes in gene_list.. if the genes are not found in the gtf_file, the fpkm value will be 0
  while line != []:
    # check if transcript or exon
    if line[2] == 'transcript':
      gene_id = ''
      gene_fpkm = ''
     
      # loop through elements on line
      for i in xrange(len(line)):
        # find gene_id
        if line[i] == 'gene_id':
          gene_id = line[i+1].strip('";')
          gene_id = gene_id.rstrip('g')
          gene_id = gene_id.rstrip('.')
          
          #print gene_id
          # check if gene_id is in gene_list
          if not gene_id in gene_list:
            gene_id = ''

        # store fpkm value
        if line[i] == 'FPKM':
          gene_fpkm = line[i+1].strip('";')
          
      # store fpkm value if gene is in gene_list
      if gene_id != '':
        try:

          if fpkm_list[gene_id] != 0:
            print 'Different Values: {} {}'.format( fpkm_list[gene_id], gene_fpkm )
          
          fpkm_list[gene_id] = gene_fpkm

        except KeyError:
          print 'Error: Key [{}] not found'.format( gene_id )

    line = gtf_fh.readline().split()

  gtf_fh.close()

  out_fh = open( outfile, 'w' )
  for i in fpkm_list:
    out_fh.write( i + '\t' + str(fpkm_list[i]) +'\n' )

  out_fh.close()


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
