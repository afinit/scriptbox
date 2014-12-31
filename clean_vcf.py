#!/usr/bin/env python
# clean vcf files by removing SNPs that are below a specified threshold and printing the remaining SNPs to a new file

import sys, os
import argparse


## Scan file to actual variant calls while printing each of the header lines to the outfile
def scan_to_snps( fh, out_fh ):
  line = fh.readline()
  while( line[0] == '#' ):
    out_fh.write( line )
    line = fh.readline()
  return line

### CLEAN VCF ###
def clean_vcf( fh, out_fh, qual, line ):
  while line != []:
    if float(line[5]) >= qual:
      out_fh.write('\t'.join(line) + '\n')
    line = fh.readline().split()

#### MAIN FUNCTION ####
def main(progname, argv):
  progname = os.path.basename( progname )
  parser = argparse.ArgumentParser( prog=progname, description='remove SNPs below a specified threshold' )
  parser.add_argument('-i, --infile', type=str, dest='infile', metavar='INFILE', required=True, help='.vcf file to be processed')
  parser.add_argument('-o, --outfile', type=str, dest='outfile', metavar='OUTFILE', help='.vcf file to be processed (if not specified, outfile is modified infile name)')
  parser.add_argument('-q, --quality', type=float, dest='qual', metavar='QUAL', required=True, help='quality threshold to process variants with')
  args = parser.parse_args(argv)

  infile = args.infile
  qual = args.qual
  if args.outfile == None:
    outfile = os.path.splitext(infile)[0]
    outfile = outfile + '.q' + str(int(qual)) + '.vcf'
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

  line = scan_to_snps( fh, out_fh )

  print 'Cleaning file'
  clean_vcf( fh, out_fh, qual, line.split() )

  fh.close()
  out_fh.close()


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
