#!/usr/bin/env python
# clean vcf files by removing SNPs that are below a specified threshold and printing the remaining SNPs to a new file

import sys, os
import argparse


### CLEAN VCF ###
def clean_file( fh, out_fh, qual ):
  line = fh.readline().split()

  while line != []:
    if float(line[3]) >= qual:
      out_fh.write('\t'.join(line) + '\n')
    line = fh.readline().split()

#### MAIN FUNCTION ####
def main(progname, argv):
  progname = os.path.basename( progname )
  parser = argparse.ArgumentParser( prog=progname, description='remove identified restriction enzymes below a specified threshold' )
  parser.add_argument('-i, --infile', type=str, dest='infile', metavar='INFILE', required=True, help='file to be processed')
  parser.add_argument('-o, --outfile', type=str, dest='outfile', metavar='OUTFILE', help='file to be written (if not specified, outfile is modified infile name)')
  parser.add_argument('-q, --quality', type=float, dest='qual', metavar='QUAL', required=True, help='quality threshold to process with')
  args = parser.parse_args(argv)

  infile = args.infile
  qual = args.qual
  if args.outfile == None:
    outfile = os.path.splitext(infile)[0]
    outfile = outfile + '.q' + str(int(qual)) + '.rec'
  else:
    outfile = args.outfile

  ## Verify input file
  if not os.path.isfile(infile):
    print 'Error: ' + infile + ' does not exist'
    sys.exit(2)

  ## Verify qual value
  if qual <= 0:
    print 'Error: quality threshold must be greater than 0'
    sys.exit()

  fh = open( infile, 'r' )
  out_fh = open( outfile, 'w' )

  print 'Cleaning file'
  clean_file( fh, out_fh, qual )

  fh.close()
  out_fh.close()


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
