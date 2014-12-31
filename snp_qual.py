#!/usr/bin/env python
# plot frequency vs quality of the snps in the vcf file given

import sys, os
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


## Scan file to actual variant calls while printing each of the header lines to the outfile
def scan_to_snps( fh ):
  line = fh.readline()
  while( line[0] == '#' ):
    line = fh.readline()
  return line

### COUNT SNPS ###
def count_snps( fh, line, freq ):
  while line != []:
    qual = int(float(line[5]))
    if len(freq) < qual:
      freq.resize((qual+1), refcheck=False)
    freq[qual] += 1
    line = fh.readline().split()

#### MAIN FUNCTION ####
def main(progname, argv):
  progname = os.path.basename( progname )
  parser = argparse.ArgumentParser( prog=progname, description='plot frequency vs SNP quality' )
  parser.add_argument('-i, --infile', type=str, dest='infile', metavar='INFILE', required=True, help='.vcf file to be processed')
  parser.add_argument('-o, --outfile', type=str, dest='outfile', metavar='OUTFILE', default='snp_qual', help='file prefix to store plot')
  parser.add_argument('-y, --ymax', type=int, dest='ymax', metavar='YMAX', default=0, help='ymax on plot')
  parser.add_argument('-x, --xmax', type=int, dest='xmax', metavar='XMAX', default=0, help='xmax on plot')
  args = parser.parse_args(argv)

  infile = args.infile
  outfile = args.outfile
  xmax = args.xmax
  ymax = args.ymax

  ## Verify input file
  if not os.path.isfile(infile):
    print 'Error: ' + infile + ' does not exist' 
    sys.exit(2)

  fh = open( infile, 'r' )

  line = scan_to_snps( fh )
  
  freq = np.zeros(40)

  count_snps( fh, line.split(), freq )

  fh.close()

  
  ###### PLOT RESULTS ######
  pref = os.path.splitext( outfile )[0]
  pdf = PdfPages( outfile + '.pdf'  )
  
  snp_plot = plt.subplot( 111 )
  plt.plot( freq )
  plt.title( 'SNP Distribution' )
  if xmax > 0:
    plt.xlim( 0, xmax )
  if ymax > 0:
    plt.ylim( 0, ymax )
  snp_plot.set_ylabel('Number of SNPs')
  snp_plot.set_xlabel('Quality')
  
  pdf.savefig()
  plt.close()
  pdf.close()


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
