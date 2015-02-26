#!/usr/bin/env python
# read a vcf file and print out average, min, and max for the following:
#   QUAL
#   MQ
#   GQ
#   DP
#   plus the frequency of genotypes
# USAGE: vcfstats.py -i infile.vcf

import os, sys
import argparse
from collections import Counter


# process stats
def proc_stats( stats, value ):
  if value < stats[2]:
    stats[2] = value
  if value > stats[3]:
    stats[3] = value
  stats[0] += value
  stats[1] += 1

################
##### MAIN #####
################
def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='read vcf file and print out min, max, and average on a few numbers',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('-i,--infile', dest='infile', metavar='INFILE', type=str, required=True, help='vcf file to provide stats for')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  infile = args.infile
  
  # check if input files exist
  if not os.path.isfile( infile ):
    print 'Error: {0} does not exist'.format( infile )
    sys.exit(1)

  # INITIALIZE
  #   index 0 = sum
  #   index 1 = count
  #   index 2 = min
  #   index 3 = max
  qual_stats = [0,0,float('inf'),-1.]
  mapq_stats = [0,0,float('inf'),-1.]
  genq_stats = [0,0,float('inf'),-1.]
  dpth_stats = [0,0,float('inf'),-1.]
  dpth_values = []
  total_snps = 0
  gt = {}

  # read in vcf entries
  print 'READING {0}'.format( infile )
  with open( infile, 'r' ) as fh:
    line = fh.readline().split()
    while line[0][0] == '#':
      line = fh.readline().split()

    while line != []:
      total_snps += 1

      # QUAL SCORE
      qual = float(line[5])
      proc_stats( qual_stats, qual )

      # parse INFO fields into a dictionary
      info_fields = dict(keypairs.split('=') for keypairs in line[7].split(';'))
  
      # MQ SCORE
      if 'MQ' in info_fields:
        proc_stats( mapq_stats, float(info_fields['MQ']) )
      else:
        mapq_stats[1] += 1
        if 0 < mapq_stats[2]:
          mapq_stats = 0

      # DEPTH
      if 'DP' in info_fields:
        proc_stats( dpth_stats, float(info_fields['DP']) )
        dpth_values.append( info_fields['DP'] )
      else:
        dpth_stats[1] += 1
        if 0 < dpth_stats[2]:
          dpth_stats = 0

      # parse genotype fields into dictionary
      try:
        gt_fields = dict(zip( line[8].split(':'), line[9].split(':') ))
      except IndexError:
        print 'Error: Invalid number of fields in vcf entry'
        sys.exit()

      # GENOTYPE
      if 'GT' in gt_fields:
        if gt_fields['GT'] in gt:
          gt[gt_fields['GT']] += 1
        else:
          gt[gt_fields['GT']] = 1
      else:
        print 'Error: GT field not present in vcf entry'
      
      # GENOTYPE QUALITY
      if 'GQ' in gt_fields:
        proc_stats( genq_stats, float(gt_fields['GQ']) )

      line = fh.readline().split()

  # calculate DP FreqDist
  print 'COUNTING DP FREQ DIST'
  dpth_freq_dist = Counter( dpth_values )
  dpth_values_sorted = sorted( dpth_freq_dist.keys(), key=lambda k: int(k) )
  with open( os.path.splitext(infile)[0] + '_dp.stats', 'w' ) as out_fh:
    for i in xrange( int(dpth_values_sorted[-1]) + 1 ):
      if str( i ) in dpth_freq_dist:
        out_fh.write( str(dpth_freq_dist[str(i)]) + '\t' )
      else:
        out_fh.write( '0' + '\t' )
    out_fh.write( '\n' )


  # print stats
  print 'VCF STATS'
  print 'SNP\tQUAL\t\t\tMQ  \t\t\tDP  \t\t\tGQ  \t\t\tGT'
  print 'total\tavg\tmin\tmax\tavg\tmin\tmax\tavg\tmin\tmax\tavg\tmin\tgt'
  
  data_out = str(total_snps) + '\t'

  # QUAL
  data_out += str(qual_stats[0]/qual_stats[1]) + '\t' + str(qual_stats[2]) + '\t' + str(qual_stats[3]) + '\t'

  # MQ SCORE
  data_out += str(mapq_stats[0]/mapq_stats[1]) + '\t' + str(mapq_stats[2]) + '\t' + str(mapq_stats[3]) + '\t'

  # DEPTH
  data_out += str(dpth_stats[0]/dpth_stats[1]) + '\t' + str(dpth_stats[2]) + '\t' + str(dpth_stats[3]) + '\t'

  # GENOTYPE QUALITY
  data_out += str(genq_stats[0]/genq_stats[1]) + '\t' + str(genq_stats[2]) + '\t'

  # GENOTYPE
  for i in sorted(gt, key=lambda k: gt[k], reverse=True):
    data_out += i + '\t' + str(gt[i]/float(total_snps)) + '\t'

  print data_out


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
