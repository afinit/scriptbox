#!/usr/bin/env python
# read in and parse a cuffdiff file into 3 files, one up, one down, and one up/down using the log2foldchange as the cutoff for those to include and only including those with yes in the significant column
# USAGE: cuffdiff_filter.py -i file.diff [-l log2foldchange] [-o output_base ]

import os, sys
import argparse


# print processed lines to file
def print_file( outfile, outlines, cuffdiff_head ):
  with open( outfile, 'w' ) as out_fh:
    out_fh.write( cuffdiff_head + '\n' )
    for line in outlines:
      out_fh.write( '\t'.join( line ) + '\n' )

################
##### MAIN #####
################
def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='process cuffdiff file into 3 smaller files',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('-i,--infile', dest='infile', metavar='INFILE', type=str, required=True, help='cuffdiff input file')
  parser.add_argument('-l,--log2foldchange', dest='foldchange', metavar='FOLDCHANGE', type=float, default=2., help='threshold value of the absolute value of log2foldchange')
  parser.add_argument('-o,--output_base', dest='outbase', metavar='OUTBASE', type=str, help='base of filenames to be output, default is infile without the extension')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  infile = args.infile
  foldchange = args.foldchange
  outbase = args.outbase

  # if outbase is not specified, create from infile
  if outbase is None:
    outbase = os.path.splitext( infile )[0]

  # check if infile exists
  if not os.path.isfile( infile ):
    print 'Error: {0} does not exist'.format( infile )
    sys.exit(1)

  # INITIALIZE
  upreg = []
  downreg = []
  filtered = []
  cuffdiff_head = ''

  print 'READING {0}'.format( infile )
  with open( infile ) as fh:
    line = fh.readline().split()
    if len(line) != 14:
      print 'Error: infile does not fit expected format'
      sys.exit()
    else:
      cuffdiff_head = '\t'.join( line )

    for line in fh:
      line = line.split()
      if line[13] != 'yes':
        continue
      elif float(line[9]) >= foldchange:
        upreg.append( line )
        filtered.append( line )
      elif float(line[9]) <= -foldchange:
        downreg.append( line )
        filtered.append( line )

  fh.close()

  print 'PRINTING OUTPUT'
  print_file( outbase + '_filt_upreg.diff', upreg, cuffdiff_head )
  print_file( outbase + '_filt_downreg.diff', downreg, cuffdiff_head )
  print_file( outbase + '_filt.diff', filtered, cuffdiff_head )

  print 'Done'

if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
