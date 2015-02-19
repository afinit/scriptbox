#!/usr/bin/env python
# sorts vcf file entries by position within their respective chromosomes
# USAGE: vcfsort.py -i infile.vcf [-o output.vcf]

import os, sys
import argparse


# replace old contig list with new contig list in header
def add_new_contigs( chrom_list, dict_list, snp_header ):
  for chrom in chrom_list:
    snp_header.append( '##contig=<ID={0},length={1}>\n'.format(chrom,dict_list[chrom]))

################
##### MAIN #####
################
def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='sort vcf file entries by position within their respective chromosomes',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('-i,--infile', dest='infile', metavar='INFILE', type=str, required=True, help='.vcf file to be processed')
  parser.add_argument('-d,--dict', dest='dictfile', metavar='DICTFILE', type=str, required=True, help='.dict file that contains proper order and length of contigs')
  parser.add_argument('-o,--outfile', dest='outfile', metavar='OUTFILE', type=str, help='output file to which results are print')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  infile = args.infile
  dictfile = args.dictfile
  outfile = args.outfile

  if outfile is None:
    outfile = os.path.splitext( infile )[0] + '_sort.vcf'

  # check if input files exist
  invalid_file = False
  if not os.path.isfile( infile ):
    print 'Error: {0} does not exist'.format( infile )
    invalid_file = True
  if not os.path.isfile( dictfile ):
    print 'Error: {0} does not exist'.format(dictfile)
    invalid_file = True
  if invalid_file:
    sys.exit(1)

  # INITIALIZE
  snp_header = []
  snp_lines = {}
  chrom_list = []
  dict_list = {}
  
  print 'READING {0}'.format( dictfile )
  with open( dictfile, 'r' ) as dict_fh:
    line = dict_fh.readline().split()

    while line != []:
      if line[0] == '@SQ':
        chrom_vals = dict(keypairs.split(':',1) for keypairs in line[1:])
        chrom_list.append(chrom_vals['SN'])
        dict_list[chrom_vals['SN']] = chrom_vals['LN']
      line = dict_fh.readline().split()

  # initialize data structure to store lines from vcf file
  snp_lines = { k: [] for k in dict_list }

  print 'READING {0}'.format( infile )
  with open( infile, 'r' ) as fh:
    line = fh.readline()

    # check 
    while line[0] == '#':
      if line[:8] == '##contig':
        add_new_contigs( chrom_list,  dict_list, snp_header )
        while line[:8] == '##contig':
          line = fh.readline()
      
      snp_header.append(line)
      line = fh.readline()

    # split line into list in prep for stage where SNPs are sorted
    line = line.split()
    snp_lines[line[0]].append(line)

    for line in fh:
      line = line.split()
      snp_lines[line[0]].append(line)

  print 'PRINTING VARIANTS'
  with open( outfile, 'w' ) as out_fh:
    for line in snp_header:
      out_fh.write( line )

    for chrom in chrom_list:
      if snp_lines[chrom] != []:
        for line in sorted( snp_lines[chrom], key=lambda k: int(k[1]) ):
          out_fh.write( '\t'.join( line ) + '\n' )

  print 'Done'


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
