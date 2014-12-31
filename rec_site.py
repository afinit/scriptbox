#!/usr/bin/env python
# usage: rec_site.py -r recognition_site_file -c contigs_file -o outfile -s vcf_file -i vcf_diff_file [-i vcf_diff_file [-i ...]]"
#   the -s vcf_file is a file used to build the genome described by those snps from the reference
#   each of the vcf_diff_file's will be used to process the contigs to discover any restriction enzymes that could cut the genome
#
# rec_sites file is formatted in a tabular file with one site per line and sites are ordered: sequence first
#  
# 
# Program Structure:
#   Read rec_sites into dict: sequence -> restriction_enzyme
#     Get max_len for rec_sites
#     Get min_len for rec_sites
#     Separate structures based on palindromic?
#   Read contigs into dict: id -> sequence
#   Open vcf_file
#   Scan to SNPs in vcf_file
#   while SNP != []
#     find string by contig_id max_len before and after pos
#     scan test each substring of lengths between min_len and max_len that include the SNP
#     ##TEST FOR MULTIPLE SNPS IN AREA???
#     
#     write matches to file

# New movement:
#   determine indels, functionality for hets(heterozygous locations)
#   snps[ chrom ][ pos ] = [[ ref, alt ]]
#   diffs[filename][ chrom ][ pos ] = [ref, alt]
#   
#   find location in snps file
#   find other snps within the max_len range
#   find all restriction_sites with snps included
#   determine number of snps within range of rec_site
#   output:   filename    scaffold    pos    restriction_enzyme
#   CONSIDER: 
#     - treatment of indels in search for rec_sites
#     - if each file has a snp but the snps are different, then any restriction enzymes found matching the reference sequence count for nothing

import signal
import sys, getopt
import os
import string

def print_usage():
  print "usage: rec_site.py -r recognition_site_file -c contigs_file -o outfile -s vcf_file -j vcf_diff_file(curr_genome) -i vcf_diff_file [-i vcf_diff_file [-i ...]]"

def signal_handler(signal, frame):
  print('SIGINT caught: Exiting')
  sys.exit(0)

# scan vcf_file to snps and return the first snp in the file
def scan_to_snps( fh ):
  line = fh.readline()
  while( line[0] == '#' ):
    line = fh.readline()
  return line.split()

# return the complement 
def complement( base ):
  if base == 'A':
    return 'T'
  elif base == 'T':
    return 'A'
  elif base == 'C':
    return 'G'
  elif base == 'G':
    return 'C'

# return reverse complement
def rev_comp( seq ):
  seq_l = list( seq )
  seq_l.reverse()
  seq_l = [ complement( x ) for x in seq_l ]
  seq = "".join( seq_l )
  return seq

def check_rec_data( seq ):
  for c in seq:
    if c not in ['A','T','C','G']:
      return False
  return True

def get_rec_sites( f ):
  rec_sites = dict()
  fh = open( f, 'r' )
  line = fh.readline().split()
  while line != []:
    rec_sites[line[1]] = line[0]
    if not check_rec_data( line[1] ):
      print "Recognition site file not in the proper format. Please make sure the file is tab delimited and has one recognition site per line with the id then sequence."
      sys.exit(2)

    # check reverse complement of current recognition site
    rev_c = rev_comp( line[1] )
    if rev_c != line[1]:
      rec_sites[ rev_c ] = line[0]
    line = fh.readline().split()
  fh.close()  
  return rec_sites
 
# get snps from vcf_file
def get_snps( fh, line ):
  # snps[ chrom ][ pos ] = [[ ref, alt ]]
  snps = {}
  # initialize current chrom in snps data structure
  snps[line[0]] = {}
  snps[line[0]][int(line[1])-1] = ( line[3], line[4] )
  chrom = line[0]

  # loop through lines of snp file
  while line != []:
    if line[0] != chrom:
      # initialize new chrom in snps data structure 
      snps[line[0]] = {}

    snps[line[0]][int(line[1])-1] = ( line[3], line[4] )
    chrom = line[0]

    # get next line
    line = fh.readline().split()

  return snps

# return contigs from reference genome in a dictionary
def get_contigs( f ):
  contigs = dict()
  id = ''
  seq = ''

  upper = str.upper
  rstrip = str.rstrip

  fh = open( f, 'r' )
  line = fh.readline()
  while line != '':
    if line[0] == '>':
      contigs[id] = seq
      id = line.split()[0][1:]
      seq = ''
    else:
      seq += upper(rstrip(line))
    line = fh.readline()

  # add last line to the last contig  
  seq += upper(rstrip(line))
  contigs[id] = seq

  # remove init item from dict
  contigs.pop( '', None )

  fh.close()
  return contigs

# read through each vcf_diff_file and return each diff in a dictionary
def get_diffs( diff_file ):
  # diffs[filename][ chrom ][ pos ] = [ref, alt]
  diffs = dict()

  # loop through each file
  for f in diff_file:
    print "\tReading: " + '\t' + f
    fh = open( f, 'r' )
    filename = os.path.basename( f )
    filename = os.path.splitext( filename )[0]
    line = scan_to_snps( fh )
    diffs[filename] = dict()

    while line != []:
      try:
        diffs[filename][line[0]].append((int(line[1])-1, float(line[5])))
      except KeyError:
        diffs[filename][line[0]] = []
        diffs[filename][line[0]].append((int(line[1])-1, float(line[5])))

      # retrieve next line in file
      line = fh.readline().split()

    fh.close()

  return diffs

# find limits for search area based
def snp_search( ref_snps, snp_pos, min_index, max_index ):
  # create list of snps within current range
  snp_list = sorted(list(set(xrange(min_index,max_index+1)).intersection(ref_snps)))
  # loop through snp_list
  for snp in snp_list:
    if snp < snp_pos and snp >= min_index:
      min_index = snp + 1
    elif snp > snp_pos and snp <= max_index:
      max_index = snp - 1
  
  # return new limits
  return min_index, max_index

# replace snp in seq
def snp_replace( seq, snp, snp_pos ):
  if seq.find( snp[0], snp_pos ) != snp_pos:
    return ''
  
  # insert snp
  seq = seq[:snp_pos] + snp[1] + seq[snp_pos+len(snp[0]):]
  return seq

# find any recognition sites containing the SNP passed
def find_sites( outstr, ref_snps, rec_sites, seq, snp, snp_pos, max_len, out_fh ):
  # replace SNP if processing genome's own diff file
  if snp != '':
    if snp[0].find(',') != -1 or snp[1].find(',') != -1:
      return
    seq = snp_replace( seq, snp, snp_pos )
    if seq == '':
      print 'Error: did not find SNP at position ' + str(snp_pos+1)
      return

  seq_len = len( seq )
  
  # find limits to account for SNPs within search area
  snp_min_index, snp_max_index = snp_search( ref_snps, snp_pos, snp_pos - (max_len-1), snp_pos + (max_len-1) )
  
  # loop over rec_sites
  for site in rec_sites:
    site_len = len(site)
    # get limits for search area
    min_index = snp_min_index if snp_pos - site_len + 1 <= snp_min_index else snp_pos - site_len + 1  # minimum index in search area
    max_index = snp_max_index if snp_pos + site_len - 1 >= snp_max_index else snp_pos + site_len - 1  # maximum index in search area

    # get search sequence from seq
    search_seq = seq[min_index:max_index+1]
    rec_site_pos = search_seq.find( site )
    if rec_site_pos != -1:
      out_fh.write( outstr + '\t' + search_seq + '\t' + str(snp_pos - min_index) + '\t' + site + '\t' + rec_sites[site] + '\n' )
      print( outstr + '\t' + search_seq + '\t' + str(snp_pos - min_index) + '\t' + str(snp) + '\t' + site + '\t' + rec_sites[site] + '\n' )

def main(argv):
  signal.signal(signal.SIGINT, signal_handler)
  snp_file = ''
  diff_file = []
  contig_file = ''
  rec_sites_file = ''
  outfile = ''
  rec_sites = dict()
  contigs = dict()
  max_rec_len = 0
  min_rec_len = 0

  # ARG PROCESSING
  try:
    opts, args = getopt.getopt(argv,"hr:c:i:j:o:s:")
  except getopt.GetoptError:
    print_usage()
    sys.exit(2)

  for opt, arg in opts:
    if opt == '-h':
      print_usage()
      sys.exit()
    elif opt == '-r':
      if not os.path.isfile( arg ):
        print "The recognition_site_file does not exist."
        print_usage()
        sys.exit(2)
      rec_sites_file = arg
    elif opt == '-s':
      if not os.path.isfile( arg ):
        print "The vcf_file does not exist."
        print_usage()
        sys.exit(2)
      snp_file = arg
    elif opt == '-i':
      if not os.path.isfile( arg ):
        print "The vcf_diff_file does not exist."
        print_usage()
        sys.exit(2)
      diff_file.append(arg)
    elif opt == '-j':
      if not os.path.isfile( arg ):
        print "The current vcf_diff_file does not exist."
        print_usage()
        sys.exit(2)
      diff_file.insert(0, arg)
    elif opt == '-c':
      if not os.path.isfile( arg ):
        print "The contigs_file does not exist."
        print_usage()
        sys.exit(2)
      contig_file = arg
    elif opt in ("-o", "--ofile"):
      outfile = arg
 
  # verify correct parameters have been selected
  if rec_sites_file == '':
    print "Recognition sites file needs to be specified"
    print_usage()
    sys.exit(2)
  if contig_file == '':
    print "Contigs file needs to be specified"
    print_usage()
    sys.exit(2)
  if outfile == '':
    print "Output file needs to be specified"
    print_usage()
    sys.exit(2)
  if snp_file == '':
    print "vcf file needs to be specified"
    print_usage()
    sys.exit(2)
  if diff_file == []:
    print "vcf_diff_file(s) need to be specified"
    print_usage()
    sys.exit(2)

  print "Reading recognition sites file"
  rec_sites = get_rec_sites( rec_sites_file )
  min_rec_len = min( [len(x) for x in rec_sites.keys()] )
  max_rec_len = max( [len(x) for x in rec_sites.keys()] )
  
  print "Reading contigs file"
  contigs = get_contigs( contig_file )
  
  # print input params
  print 'SNP file is ' + snp_file
  print 'Output file is ' + outfile
  for f in diff_file:
    print 'diff_file: ' + f

  # declare variables
  in_fh = open( snp_file, 'r' ) # file handler for input file
  file_line = scan_to_snps( in_fh ) # current line
  
  print "Reading SNP file for current genome"
  ref_snps = get_snps( in_fh, file_line )
  in_fh.close()

  print "Reading SNP difference files"
  diffs = get_diffs( diff_file )
  out_fh = open( outfile, 'w' ) # file handler for output file
  rec_sites_keys = rec_sites.keys()

  native_snpfile = os.path.basename( diff_file[0] )
  native_snpfile = os.path.splitext( native_snpfile )[0]

  for filename in diffs:
    for chrom in sorted( diffs[filename] ):
      for (pos, qual) in sorted( diffs[filename][chrom] ):
        print "Processing: filename: " + filename + " chrom: " + chrom + " pos: " + str(pos)
        ### snp_seq = diffs[chrom][pos]
        if( filename == native_snpfile ):
          find_sites( filename + '\t' + chrom + '\t' + str(pos+1) + '\t' + str(qual), ref_snps[chrom], rec_sites, contigs[chrom], ref_snps[chrom][pos], pos, max_rec_len, out_fh )
        else:
          find_sites( filename + '\t' + chrom + '\t' + str(pos+1) + '\t' + str(qual), ref_snps[chrom], rec_sites, contigs[chrom], '', pos, max_rec_len, out_fh )

  # close files
  out_fh.close()

if __name__ == "__main__":
  main(sys.argv[1:])
