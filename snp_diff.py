#!/usr/bin/env python
#usage: snp_diff.py -i infile -o outtag
#

import sys, getopt
import os

def scan_to_snps( fh, out_fh ):
  chrom = []
  line = fh.readline()
  while( line[0] == '#' ):
    out_fh.write( line )
    if line[:8] == '##contig':
      id = line[10:-1]
      id = dict(x.split("=") for x in id.split(","))
      chrom.append( id['ID'] )
    line = fh.readline()
  return [line, chrom]

def pull_current_chrom( fh, line, chrom ):
  snps = dict()
  while line != [] and line[0] == chrom:
    snps[line[1]] = line
    line = fh.readline().split()
  return [line, snps]

def main(argv):
  infile = []
  outtag = ''

  try:
    opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
  except getopt.GetoptError:
    print 'snp_diff.py -i <infile> -i <infile> [-i <infile>...] -o <outtag>'
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print 'snp_diff.py -i <infile> -i <infile> [-i <infile>...] -o <outtag>'
      sys.exit()
    elif opt in ("-i", "--ifile"):
      infile.append(arg)
    elif opt in ("-o", "--ofile"):
      outtag = '.' + arg

  # print input params
  for f in infile:
    print 'Input file is ', f

  # declare variables
  files = [] # file handlers for input files
  out_fh = [] # file handlers for output files
  file_line = [] # current line for each file
  files_left = [] # indices of the files on the current chrom
  curr_chrom = ''
  count_diff = []
  total = []
  count_same = 0
  chrom = []

  #initialize variables
  for f in infile:
    files.append( open( f, "r" ) )
    out_fh.append( open( os.path.splitext(f)[0] +  outtag + '.diff.vcf', 'w' ) )
    count_diff.append(0)
    total.append(0)

  for i in range(len(files)):
    file_scan = scan_to_snps( files[i], out_fh[i] )
    chrom = file_scan[1]
    file_line.append( file_scan[0].split() )

  # process each chrom in the file
  for curr_chrom in chrom:
    print "Processing: " + curr_chrom
    # set files_left variable and curr_pos variable for iterating
    snps = []
    snps_keys = []
    key_list = []

    # pull data for current chrom
    for i in xrange(len(files)):
      pull_data = pull_current_chrom( files[i], file_line[i], curr_chrom )
      file_line[i] = pull_data[0]
      snps.append( pull_data[1] )
      snps_keys.append( pull_data[1].keys() )

    # create list of all represented postions
    for snp_dict in snps:
      key_list = list( set( snp_dict ) - set( key_list ) ) + key_list

    key_list = sorted( key_list )

    # loop through each represented position
    for curr_pos in key_list:
      
      # for tracking indices matching the current position
      pos_match = []

      # loop through each file left on the current chrom
      for i in xrange(len(snps)):
        try:
          snps[i][curr_pos]
          pos_match.append(i)
        except KeyError:
          pass

      # check if all match at this point
      # if true, increase count_same, and total for each file and check that each snp is the same
      # check that each snp is the same
      if len( pos_match ) == len( files ):
        diff_flag = False
        alt_bp = snps[0][curr_pos][4]
        for f in snps:
          if f[curr_pos][4] != alt_bp:
            count_diff = [x+1 for x in count_diff]
            for i in range(len(snps)):
              out_fh[i].write( '\t'.join( snps[i][curr_pos] ) + '\n' )

            diff_flag = True
            break

        total = [x+1 for x in total]
        if not diff_flag:
          count_same += 1
      
      # if false, increase total and count_diff for each file included in the list
      else:
        for i in pos_match:
          total[i] += 1
          count_diff[i] += 1
          out_fh[i].write( '\t'.join( snps[i][curr_pos] ) + '\n' )

  for i in range(len(files)):
    if total[i] > 0:
      print infile[i] + " diff: " + str( count_diff[i]) + " (" + str(float(count_diff[i])/total[i]*100) + "%)"
    print infile[i] + " total: " + str( total[i] )

  print "same: " + str(count_same)

  # close files
  for f in files:
    f.close()
  
  for fh in out_fh:
    fh.close()

if __name__ == "__main__":
  main(sys.argv[1:])
