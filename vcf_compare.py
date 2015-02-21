#!/usr/bin/env python
# reads in and stores multiple vcf files
# USAGE: vcf_compare.py -r ref.fa INFILES.vcf

import os, sys
import argparse


################
##### MAIN #####
################
def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='reads in and stores multiple vcf files',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('infile', metavar='INFILE', type=str, nargs='+', help='vcf files to be processed')
  parser.add_argument('-r','--ref', dest='ref_file', metavar='REF_FILE', type=str, required=True, help='reference file which SNPs were called against, must be indexed with Picard Tool\'s CreateSequenceDictionary')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  infile = args.infile
  ref_file = args.ref_file

  if infile is None or len(infile) < 2:
    print 'Error: Must supply at least two vcf files to compare'
    sys.exit(1)
  
  # check if input files exist
  invalid_file = False
  if not os.path.isfile( gff_file ):
    print 'Error: {0} does not exist'.format( gff_file )
    invalid_file = True
  if not os.path.isfile( vcf_file ):
    print 'Error: {0} does not exist'.format( vcf_file )
    invalid_file = True
  if invalid_file:
    sys.exit(1)

  # set outfile if not specified
  if outfile == '':
    outfile = os.path.splitext( vcf_file )[0] + '_lcfltrd.vcf'

  # INITIALIZE
  lc_ranges = {}
  variants = {}
  variants_proc = []
  vcf_header = []

  # read in lc_ranges from gff file
  print 'READING {0}'.format( gff_file )
  with open( gff_file, 'r' ) as gff_fh:
    line = gff_fh.readline().split()
    while line[0][0] == '#':
      line = gff_fh.readline().split()

    while line != []:
      if not line[0] in lc_ranges:
        lc_ranges[line[0]] = []

      lc_ranges[line[0]].append( (int(line[3]), int(line[4])+1) )
      line = gff_fh.readline().split()

  # read in variants from vcf file
  print 'READING {0}'.format( vcf_file )
  with open( vcf_file, 'r' ) as vcf_fh:
    line = vcf_fh.readline().split()
    while line[0][0] == '#':
      vcf_header.append( line )
      line = vcf_fh.readline().split()
    
    while line != []:
      if not line[0] in variants:
        variants[line[0]] = []

      variants[line[0]].append( line )
      line = vcf_fh.readline().split()

  # loop through variants
  print 'PROCESSING VARIANTS'
  for chrom in variants:
    if chrom in lc_ranges:
      process_lowcmplx( lc_ranges[chrom], variants[chrom], variants_proc )
      pass
    else:
      variants_proc += variants[chrom]

  print 'WRITING OUTPUT TO FILE'
  with open( outfile, 'w' ) as out_fh:
    for line in vcf_header:
      out_fh.write( '\t'.join( line ) + '\n' )

    for line in variants_proc:
      out_fh.write( '\t'.join( line ) + '\n' )

  print 'Done'


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
