#!/usr/bin/env python
# reads in and stores multiple vcf files
# USAGE: vcf_compare.py -r ref.fa INFILES.vcf

import os, sys
import argparse


# read dict file
def read_dict_file(dict_file, dict_list, chrom_list):
    print 'READING {0}'.format( dict_file )
    with open( dict_file, 'r' ) as dict_fh:
        line = dict_fh.readline().split()
    
        while line != []: 
        if line[0] == '@SQ':
            chrom_vals = dict(keypairs.split(':',1) for keypairs in line[1:])
            chrom_list.append(chrom_vals['SN'])
            dict_list[chrom_vals['SN']] = chrom_vals['LN']
        line = dict_fh.readline().split()
    print 'done reading'

# read in variants from vcf file
def read_vcf_file(vcf_file, SNPs, file_index):
    print 'READING {0}'.format( vcf_file )
    with open( vcf_file, 'r' ) as vcf_fh:
        line = vcf_fh.readline()
    
        # loop through header lines
        while line[0] == '#':
            line = vcf_fh.readline()
    
            line = line.split()

        # loop through variants
        while line != []:
            SNPs[int(line[1])][file_index] = line[4]
            line = vcf_fh.readline().split()

# read through file to find offsets for each chrom
def get_offsets( chrom_offset, vcf_file, chrom_list ):
    if not os.path.isfile( vcf_file + '.chromi' ) or os.path.getmtime( vcf_file ) < os.path.getmtime( vcf_file + '.chromi' ):
        print 'BUILDING OFFSET CHROM FOR {0}'.format( vcf_file )
        chrom_offset = { k: None for k in chrom_list }
        
        # Read in the file once and build a list of line offsets
        offset = 0
        current_chrom = ''
        for line in fh:
            if line[0] != '#' and current_chrom != line.split()[0]:
                chrom_offset[current_chrom] = offset
                current_chrom = line.split()[0]
            offset += len(line)
    
        # write offsets index for each chrom in the vcf file
        with open( vcf_file + '.chromi', 'w') as out_fh:
          out_fh.write( str(chrom_offset) )
    else:
        with open( vcf_file + '.chromi', 'r' ) as fh:
            chrom_offset = {}
            chrom_offset = eval(fh.readline())

# build SNPs data structure and save it to file
def build_snp_structure( chrom_list ):
    for chrom in chrom_list:
        




    #for chrom in chrom_list:


    #  print 'creating SNPs data structure'
    # initialize data structure to store lines from vcf file
    #  SNPs = { p:  [None for x in xrange(len(infile))]  for p in xrange(int(dict_list[chrom]))}



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
    dict_file = os.path.splitext(ref_file)[0] + '.dict'

    if infile is None or len(infile) < 2:
        print 'Error: Must supply at least two vcf files to compare'
        sys.exit(1)
  
    # check if input files exist
    invalid_file = False
    for f in infile:
        if not os.path.isfile( f ):
            print 'Error: {0} does not exist'.format( f )
            invalid_file = True
        if not os.path.isfile( ref_file ):
            print 'Error: {0} does not exist'.format( ref_file )
            invalid_file = True
        if not os.path.isfile( dict_file ):
            print 'Error: {0} does not exist. The reference file (ref.fa) must have a corresponding .dict file (ref.dict) created by Picard Tool\'s CreateSequenceDictionary'.format( ref_file )
            invalid_file = True
        if invalid_file:
            sys.exit(1)

    # INITIALIZE
    SNPs = {}
    chrom_list = []
    dict_list = {}
    chrom_offset = {}
  
    # read .dict file
    read_dict_file(dict_file, dict_list, chrom_list)

    for vcf_file in infile:
        print vcf_file
    
        get_offsets( chrom_offset, vcf_file, chrom_list )
        with open( vcf_file, 'r' ) as vcf_fh:

    #  print 'done creating SNPs data structure'

    # setup file index numbers
    file_indnums = dict( zip( infile, xrange( len( infile ) ) ))
    
    # loop through vcf files
    #  for vcf_file in infile:
    #    read_vcf_file(vcf_file, SNPs, file_indnums[vcf_file])


    print 'Done'


if __name__=='__main__':
    main(sys.argv[0], sys.argv[1:])
