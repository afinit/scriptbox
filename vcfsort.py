#!/usr/bin/env python
# sorts vcf file entries by position within their respective chromosomes
# USAGE: vcfsort.py -i infile.vcf [-o output.vcf]

import os, sys
import argparse


################
##### MAIN #####
################
def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='sort vcf file entries by position within their respective chromosomes',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('-i,--infile', dest='infile', metavar='INFILE', type=str, required=True, help='vcf file to be processed')
  parser.add_argument('-o,--outfile', dest='outfile', metavar='OUTFILE', type=str, help='output file to which results are print')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  infile = args.infile
  outfile = args.outfile

  if outfile is None:
    outfile = os.path.splitext( infile )[0] + '_sort.vcf'

  # check if input files exist
  if not os.path.isfile( infile ):
    print 'Error: {0} does not exist'.format( infile )
    sys.exit(1)

  # INITIALIZE
  snp_header = []
  snp_lines = {}
  chrom_list = []

  print 'READING {0}'.format( infile )
  with open( infile, 'r' ) as fh:
    line = fh.readline()

    while line[0] == '#':
      snp_header.append(line)
      line = fh.readline()

    line = line.split()
    while line != []:
      if not line[0] in snp_lines:
        snp_lines[line[0]] = []
        chrom_list.append(line[0])

      snp_lines[line[0]].append( line )
      
      line = fh.readline().split()


  print 'PRINTING VARIANTS'
  with open( outfile, 'w' ) as out_fh:
    for line in snp_header:
      out_fh.write( line )

    for chrom in chrom_list:
      for line in sorted( snp_lines[chrom], key=lambda k: int(k[1]) ):
        out_fh.write( '\t'.join( line ) + '\n' )

  print 'Done'


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
