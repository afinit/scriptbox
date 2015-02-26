#!/usr/bin/env python
# accepts a gff annotation file that was output from RepeatMasker (any gff would work, but this is its designed purpose) and a vcf file, it will then remove any SNPs that fall within the ranges specified in the gtf
#   file and print the SNPs back to a different vcf file unless the original file is specified in the output option
# USAGE: vcf_lowcmplx_fltr.py -v infile.vcf -g infile.gff [-o output.vcf]

import os, sys
import argparse


# remove SNPs in low complexity regions 
def process_lowcmplx( lc_ranges, variants, variants_proc ):
  lc_index = 0
  lc_end = False

  for v in variants:
    if not lc_end:
      if int(v[1]) < lc_ranges[lc_index][0]:
        variants_proc.append( v )
      elif int(v[1]) > lc_ranges[lc_index][1]:
        variants_proc.append( v )
        lc_index += 1
        if lc_index == len( lc_ranges ):
          lc_end = True
        pass
    else:
      variants_proc.append( v )

################
##### MAIN #####
################
def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='process out SNPs from vcf located in ranges specified by gff file',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('-g,--gff_file', dest='gff_file', metavar='GFF_FILE', type=str, required=True, help='annotation file containing the list of ranges specified as low complexity by RepeatMasker')
  parser.add_argument('-v,--vcf_file', dest='vcf_file', metavar='VCF_FILE', type=str, required=True, help='file containing SNPs and other variants to be processed')
  parser.add_argument('-o,--outfile', dest='outfile', metavar='OUTFILE', type=str, default='', help='output file to write to')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  gff_file = args.gff_file
  vcf_file = args.vcf_file
  outfile = args.outfile
  
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
