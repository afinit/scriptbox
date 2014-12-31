#!/usr/bin/env python
#usage: snp_dist.py -i infile -o outfile [-q cutoff] [-s snp_cutoff]
#

import sys, getopt
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def print_usage():
  print 'snp_dist.py -i <infile> [-i <infile>] -o <outfile> [-q quality_cutoff]'
  sys.exit(1)

def scan_to_snps( fh ):
  chrom = []
  length = {}
  line = fh.readline()
  while( line[0] == '#' ):
    if line[:8] == '##contig':
      id = line[10:-2]
      id = dict(x.split("=") for x in id.split(","))
      chrom.append( id['ID'] )
      length[id['ID']] = id['length']
    line = fh.readline()
  return [line, chrom, length]

def pull_current_chrom( fh, line, chrom, snp_pos, cutoff ):
  #print "chrom: " + chrom + " line: " + str(line)
  while line != [] and line[0] == chrom:
    if float(line[5]) >= cutoff:
      snp_pos.append(int(line[1]))
    line = fh.readline().split()
  #print "chrom: " + chrom + " snps: " + str(len(snps))
  snp_pos = list(set(snp_pos))
  return [line, snp_pos]

def main(argv):
  infile = []
  outfile = 'snp_dist.out'
  cutoff = 0.0
  snp_cutoff = 0
  process_jump = 500
  process_length = 1000

  try:
    opts, args = getopt.getopt(argv,"hq:s:j:l:i:o:",["ifile=","ofile="])
  except getopt.GetoptError:
    print_usage()

  for opt, arg in opts:
    if opt == '-h':
      print_usage()
    elif opt in ("-q"):
      cutoff = float(arg)
    elif opt in ("-s"):
      snp_cutoff = int(arg)
    elif opt in ("-j"):
      process_jump = int(arg)
    elif opt in ("-l"):
      process_length = int(arg)
    elif opt in ("-i", "--ifile"):
      infile.append(arg)
    elif opt in ("-o", "--ofile"):
      outfile = arg
     
  # print input params
  print 'Input file is ', infile
  print 'Quality cutoff is ', cutoff

  # declare variables
  fh = [] # file handler
  out_fh = '' # file handler for output file
  file_line = [] # current line for file
  files_left = [] # indices of the files on the current chrom
  chrom_length = {}
  curr_chrom = ''
  num_snps = 0
  bases = 0
  chrom = []
  snps = []
  snps_keys = {}
  snp_pos = {}
  region_count = 0
  percentiles = [0]*20
  snpdist = []

  #initialize variables
  for f in infile:
    fh.append( open( f, "r" ) )

  out_fh = open( outfile, 'w' )
  out_fh.write( 'chrom\tseq_len\tpos\tnum_snps ### jump: ' + str( process_jump ) + '\n' )

  # scan to snps
  for i in xrange( len(fh) ):
    file_scan = scan_to_snps( fh[i] )
    chrom = file_scan[1]
    file_line.append( file_scan[0].split() )
    snps.append({})
  chrom_length = file_scan[2]

  # process each chrom in the file
  for curr_chrom in chrom:
    #print "Processing: " + curr_chrom + ' length: ' + chrom_length[curr_chrom]
    snp_pos[curr_chrom] = []
    
    # get data for current current chrom from each file
    for i in xrange( len(fh) ):
      # pull data
      pull_data = pull_current_chrom( fh[i], file_line[i], curr_chrom, snp_pos[curr_chrom], cutoff )
      file_line[i] = pull_data[0]
      snp_pos[curr_chrom] = pull_data[1]

    # uniquify snp_pos list
    snp_pos[curr_chrom] = set( snp_pos[curr_chrom] )

    curr_len = int(chrom_length[curr_chrom])
    # process data 
    for i in xrange( 0, curr_len, process_jump ):
      stop_pos = i + process_jump
      if stop_pos > curr_len:
        stop_pos = curr_len

      local_snps = len(snp_pos[curr_chrom].intersection(xrange(i,stop_pos)))
      
      percentiles[int(float(local_snps)/process_length*200)] += 1

      if local_snps >= snp_cutoff:
        region_count += 1
        out_fh.write( curr_chrom + '\t' + str( process_length ) + '\t' + str(i) + '\t' + str(local_snps) + '\n' )
        snpdist.append( local_snps )
      else:
        snpdist.append( 0 )

  # close files
  for i in range( len(fh) ):
    fh[i].close()
  out_fh.close()
  print 'Regions counted: ' + str(region_count)

  for i in xrange(20):
    print str(float(i*5)/10) + '%: ' + str(percentiles[i])

  ###### PLOT RESULTS ######
  pref = os.path.splitext( outfile )[0]
  pdf = PdfPages( pref+".plot.pdf" )
  
  snp_plot = plt.subplot( 111 )
  plt.plot( snpdist, 'ro' )
  plt.title( 'SNP Distribution' )
  snp_plot.set_ylabel('Number of SNPs')
  snp_plot.set_xlabel('Genome Regions (scaffolds concatenated in ' + str(process_length) + 'bp bins)')
  
  pdf.savefig()
  plt.close()
  pdf.close()
  
if __name__ == "__main__":
  main(sys.argv[1:])
