#!/usr/bin/env python
# provides several filter methods:
#   MQ low limit
#   DP low limit
#   low complexity specified by gff file created by RepeatMasker (any gff would work, but this is its designed purpose) and a vcf file, it will then remove any SNPs that fall within the ranges specified in the gtf
#     file and print the SNPs back to a different vcf file unless the original file is specified in the output option
#   GQ low limit
# the remaining SNPs are then printed to the output file specified by -o
# USAGE: vcf_filter.py -i infile.vcf [-m MQ_limit] [-d DP_limit] [-q GQ_limit] [-g infile.gff] [-t genotype] -o output.vcf

import os, sys
import argparse


# retrieve variants from vcf file 
def get_variants( vcf_file, vcf_header, variants, chrom_list ):
    # read in variants from vcf file
    print 'READING {0}'.format( vcf_file )
    with open( vcf_file, 'r' ) as vcf_fh:
        line = vcf_fh.readline()
        while line[0] == '#':
            vcf_header.append( line )
            if line[:8] == '##contig':
                n_start = line.find( 'ID=' ) + 3
                n_end = line.find( ',', n_start )
                if n_end == -1:
                    n_end = line.find( '>', n_start )
                chrom_list.append( line[n_start:n_end] )
            line = vcf_fh.readline()
    
        line = line.split()
        while line != []:
            if not line[0] in variants:
                variants[line[0]] = []

            variants[line[0]].append( line )
            line = vcf_fh.readline().split()


# retrieve low complexity ranges from gff file
def get_lc_ranges( gff_file, lc_ranges ): 
    # read in lc_ranges from gff file
    print 'READING {0}'.format( gff_file )
    with open( gff_file, 'r' ) as gff_fh:
        line = gff_fh.readline().split()
        while line[0][0] == '#':
            line = gff_fh.readline().split()

        while line != []:
            if not line[0] in lc_ranges:
                lc_ranges[line[0]] = set()

            lc_ranges[line[0]] = lc_ranges[line[0]].union( xrange(int(line[3]), int(line[4])+1) )
            line = gff_fh.readline().split()

# process low complexity regions
def process_lowcomplex( variants, lc_ranges, variants_proc ):
    print 'PROCESSING LOW COMPLEXITY REGIONS'
    for chrom in variants:
        variants_proc[chrom] = []

        # if chrom has a low complexity range, check all of the snps against it/them
        if chrom in lc_ranges:
            for v in variants[chrom]:
                if not int(v[1]) in lc_ranges[chrom]:
                    variants_proc[chrom].append( v )
        # if chrom does not have a low complexity range, add all of the snps from chrom
        else:
            variants_proc[chrom] = variants[chrom]

# filter out SNPs not meeting criteria of data points contained in INFO field
def process_info_field( variants, variants_proc, mq, dp ):
    print 'PROCESSING INFO FIELD CRITERIA'
    
    # set values that weren't specified to neutral values that won't filter anything out
    if mq is None:
        mq = 0
    elif dp is None:
        dp = 0

    # initialize KeyError counter
    keyerr = 0

    # loop through chroms
    for chrom in variants:
        variants_proc[chrom] = []
        
        # loop through variants of current chrom
        for v in variants[chrom]:
            # parse INFO field
            info_fields = dict(keypairs.split('=') for keypairs in v[7].split(';'))

            # guard against keys not existing in INFO field.. altho this should never happen
            try:
                if float(info_fields['MQ']) >= mq and float(info_fields['DP']) >= dp:
                    variants_proc[chrom].append( v )
            except KeyError:
                keyerr += 1

    # alert user to keyerror occurrence 
    if keyerr:
        print 'Warning: {0} records were missing at least one key in the INFO field data'.format(keyerr)


# write filtered SNPs to file
def write_snps( outfile, vcf_header, variants, chrom_list ):
    print 'WRITING OUTPUT TO FILE'
    with open( outfile, 'w' ) as out_fh:
        for line in vcf_header:
            out_fh.write( line )

        for chrom in chrom_list:
            for line in variants[chrom]:
                out_fh.write( '\t'.join( line ) + '\n' )


# remove SNPs failing to meet the limits set at runtime 
def run_filters( variants, gff_file, mq, dp, gq, gt ):
    if gff_file:
        # initialize holding variables for low complexity processing
        lc_ranges = {}
        variants_proc = {}

        get_lc_ranges( gff_file, lc_ranges )
        process_lowcomplex( variants, lc_ranges, variants_proc )
        variants = variants_proc

    # run filters on data in INFO field
    if mq or dp:
        variants_proc = {}
        process_info_field( variants, variants_proc, mq, dp )
        variants = variants_proc

    if gq:
        pass
    
    if gt:
        pass

    return variants


################
##### MAIN #####
################
def main( prog_name, argv ):
    # ARG PROCESSING
    parser = argparse.ArgumentParser( prog=prog_name, description='process out SNPs from vcf located in ranges specified by gff file, only filters specified will be processed',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument('-i','--vcf_file', dest='vcf_file', metavar='VCF_FILE', type=str, required=True, help='file containing SNPs and other variants to be processed')
    parser.add_argument('-g','--gff_file', dest='gff_file', metavar='GFF_FILE', type=str, help='annotation file containing the list of ranges specified as low complexity by RepeatMasker')
    parser.add_argument('-m','--mq', dest='mq', metavar='MQ', type=float, help='minimum mapping quality filter')
    parser.add_argument('-d','--dp', dest='dp', metavar='DP', type=float, help='minimum DP value filter')
    parser.add_argument('-q','--gq', dest='gq', metavar='GQ', type=float, help='minimum GQ value filter')
    parser.add_argument('-t','--gt', dest='gt', metavar='GT', type=str, help='genotype to filter out in gt format( eg, 0/1 )')
    parser.add_argument('-o','--outfile', dest='outfile', metavar='OUTFILE', type=str, help='output file to write to')
    args = parser.parse_args(argv)

    # SET VARIBLES FROM COMMAND LINE ARGS
    vcf_file = args.vcf_file
    gff_file = args.gff_file
    mq = args.mq
    dp = args.dp
    gq = args.gq
    gt = args.gt
    outfile = args.outfile
  
    # check if input files exist
    invalid_file = False
    if gff_file and not os.path.isfile( gff_file ):
        print 'Error: {0} does not exist'.format( gff_file )
        invalid_file = True
    if not os.path.isfile( vcf_file ):
        print 'Error: {0} does not exist'.format( vcf_file )
        invalid_file = True
    if invalid_file:
        sys.exit(1)

    # set outfile if not specified
    if outfile is None:
        outfile = os.path.splitext( vcf_file )[0] + '_filter.vcf'

    # INITIALIZE
    variants = {}
    vcf_header = []
    chrom_list = []

    # READ VARIANTS
    get_variants( vcf_file, vcf_header, variants, chrom_list )

    # PROCESS FILTERS
    variants = run_filters( variants, gff_file, mq, dp, gq, gt )

    # WRITE FILTERED VARIANTS
    write_snps( outfile, vcf_header, variants, chrom_list )


    print 'Done'


if __name__=='__main__':
    main(sys.argv[0], sys.argv[1:])
