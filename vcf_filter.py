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


##################
## FILTER CLASS ##
##################
class Filter:
    def __init__(self, gff_file, indel, mq, dp, gq, gt, no_count ):
        # set initialize user set variables
        self.gff_file = gff_file
        self.indel = indel
        self.mq = mq
        self.dp = dp
        self.gq = gq
        self.gt = gt

        # variant file variables
        self.variants = {}
        self.vcf_header = []
        self.chrom_list = []

        # low complexity ranges from 
        self.lc_ranges = None

        # create variables for counting removed SNPs
        self.lc_rem_count = None
        self.indel_rem_count = None
        self.mq_rem_count = None
        self.dp_rem_count = None
        self.gq_rem_count = None
        self.gt_rem_count = None

        # set variables to count
        if not no_count:
            # low complexity removal counting
            if gff_file:
                self.lc_rem_count = 0

            # indel removal counting
            if indel:
                self.indel_rem_count = 0

            # MQ removal counting
            if mq:
                self.mq_rem_count = 0

            # DP removal counting
            if dp:
                self.dp_rem_count = 0

            # GQ removal counting
            if gq:
                self.gq_rem_count = 0

            # GT removal counting
            if gt:
                self.gt_rem_count = 0


    ################
    ## FILE READING

    # retrieve variants from vcf file 
    def get_variants( self, vcf_file ):
        # read in variants from vcf file
        print 'READING {0}'.format( vcf_file )
        with open( vcf_file, 'r' ) as vcf_fh:
            line = vcf_fh.readline()
            while line[0] == '#':
                self.vcf_header.append( line )
                if line[:8] == '##contig':
                    n_start = line.find( 'ID=' ) + 3
                    n_end = line.find( ',', n_start )
                    if n_end == -1:
                        n_end = line.find( '>', n_start )
                    self.chrom_list.append( line[n_start:n_end] )
                line = vcf_fh.readline()
        
            line = line.split()
            while line != []:
                if not line[0] in self.variants:
                    self.variants[line[0]] = []

                self.variants[line[0]].append( line )
                line = vcf_fh.readline().split()


    # retrieve low complexity ranges from gff file
    def get_lc_ranges( self ): 
        # read in lc_ranges from gff file
        print 'READING {0}'.format( self.gff_file )
        with open( self.gff_file, 'r' ) as gff_fh:
            line = gff_fh.readline().split()
            while line[0][0] == '#':
                line = gff_fh.readline().split()

            while line != []:
                if not line[0] in self.lc_ranges:
                    self.lc_ranges[line[0]] = set()

                self.lc_ranges[line[0]] = self.lc_ranges[line[0]].union( xrange(int(line[3]), int(line[4])+1) )
                line = gff_fh.readline().split()

    #####################
    ## FILTER PROCESSING 

    # process low complexity regions
    def process_lowcomplex( self, variants_proc ):
        print 'PROCESSING LOW COMPLEXITY REGIONS'
        for chrom in self.variants:
            variants_proc[chrom] = []

            # if chrom has a low complexity range(s), check all of the snps against it/them
            if chrom in self.lc_ranges:
                for v in self.variants[chrom]:
                    if not int(v[1]) in self.lc_ranges[chrom]:
                        variants_proc[chrom].append( v )
                    elif not self.lc_rem_count is None:
                        self.lc_rem_count += 1
            # if chrom does not have a low complexity range, add all of the snps from chrom
            else:
                variants_proc[chrom] = self.variants[chrom]

    # process indel filter (remove indels)
    def process_indel( self, variants_proc ):
      print 'REMOVING INDELS'
      for chrom in self.variants:
          variants_proc[chrom] = []

          for v in self.variants[chrom]:
              if len(v[3]) == 1 and len(v[4]) == 1:
                  variants_proc[chrom].append( v )
              else:
                  if not self.indel_rem_count is None:
                      self.indel_rem_count += 1

    # filter out SNPs not meeting criteria of data points contained in INFO field
    def process_info_field( self, variants_proc ):
        print 'PROCESSING INFO FIELD CRITERIA'
        
        # set values that weren't specified to neutral values that won't filter anything out
        if self.mq is None:
            self.mq = 0
        elif self.dp is None:
            self.dp = 0

        # initialize KeyError counter
        keyerr = 0

        # loop through chroms
        for chrom in self.variants:
            variants_proc[chrom] = []
            
            # loop through variants of current chrom
            for v in self.variants[chrom]:
                # parse INFO field
                info_fields = dict(keypairs.split('=') for keypairs in v[7].split(';'))

                # guard against keys not existing in INFO field.. altho this should never happen
                try:
                    if float(info_fields['MQ']) >= self.mq and float(info_fields['DP']) >= self.dp:
                        variants_proc[chrom].append( v )
                    else:
                        # increase count on appropriate count variable since this variant was not included
                        if not self.mq_rem_count is None and float(info_fields['MQ']) < self.mq:
                            self.mq_rem_count += 1
                        elif not self.dp_rem_count is None:
                            self.dp_rem_count += 1
                except KeyError:
                    # increase count on appropriate count variable since this variant was not included
                    if not self.mq_rem_count is None:
                        self.mq_rem_count += 1
                    elif not self.dp_rem_count is None:
                        self.dp_rem_count += 1
                    keyerr += 1

        # alert user to keyerror occurrence 
        if keyerr:
            print 'Warning: {0} records were missing at least one key in the INFO field data'.format(keyerr)

    # filter out SNPs not meeting the criteria of data points contained in genotype information
    def process_genotype_field( self, variants_proc ):
        print 'PROCESSING GENOTYPE CRITERIA'

        # set values that aren't specified to neutral values
        if self.gq is None:
            self.gq = 0
        elif self.gt is None:
            self.gt = []

        # initialize KeyError counter
        keyerr = 0

        # loop through chroms
        for chrom in self.variants:
            variants_proc[chrom] = []

            # loop through variants of current chrom
            for v in self.variants[chrom]:
                # parse GT field
                gt_fields = dict(zip( v[8].split(':'), v[9].split(':') ))

                # guard against keys not existing in INFO field.. altho this should never happen
                try:
                    if float(gt_fields['GQ']) >= self.gq and not gt_fields['GT'] in self.gt:
                        variants_proc[chrom].append(v)
                    else:
                        # increase count on appropriate count variable since this variant was not included
                        if not self.gq_rem_count is None and float(gt_fields['GQ']) < self.gq:
                            self.gq_rem_count += 1
                        elif not self.gt_rem_count is None:
                            self.gt_rem_count += 1
                except KeyError:
                    # increase count on appropriate count variable since this variant was not included
                    if not self.gq_rem_count is None:
                        self.gq_rem_count += 1
                    elif not self.gt_rem_count is None:
                        self.gt_rem_count += 1
                    keyerr += 1

        # alert user to keyerror occurrence 
        if keyerr:
            print 'Warning: {0} records were missing at least one key in the Genotype field data'.format(keyerr)

    ##########
    ## OUTPUT 

    # write filtered SNPs to file
    def write_snps( self, outfile ):
        print 'WRITING OUTPUT TO FILE'
        with open( outfile, 'w' ) as out_fh:
            for line in self.vcf_header:
                out_fh.write( line )

            for chrom in self.chrom_list:
                for line in self.variants[chrom]:
                    out_fh.write( '\t'.join( line ) + '\n' )

    # print removal counts to the screen
    def write_rem_counts( self ):
        print 'VARIANTS REMOVED (These counts are made in order of filter application, so no variants are counted 2x):'
        
        total_removed = 0

        if not self.lc_rem_count is None:
            total_removed += self.lc_rem_count
            print '\tLow Complexity: {0}'.format( self.lc_rem_count )
        if not self.indel_rem_count is None:
            total_removed += self.indel_rem_count
            print '\tIndels: {0}'.format( self.indel_rem_count )
        if not self.mq_rem_count is None:
            total_removed += self.mq_rem_count
            print '\tMapping Quality: {0}'.format( self.mq_rem_count )
        if not self.dp_rem_count is None:
            total_removed += self.dp_rem_count
            print '\tDepth: {0}'.format( self.dp_rem_count )
        if not self.gq_rem_count is None:
            total_removed += self.gq_rem_count
            print '\tGenotype Quality: {0}'.format( self.gq_rem_count )
        if not self.gt_rem_count is None:
            total_removed += self.gt_rem_count
            print '\tGenotype: {0}'.format( self.gt_rem_count )

        print '\tTotal Removed: {0}'.format( total_removed )

    ##################
    ## MANAGE FILTERS

    # remove SNPs failing to meet the limits set at runtime 
    def run_filters( self ):
        if self.gff_file:
            # initialize holding variables for low complexity processing
            self.lc_ranges = {}
            variants_proc = {}

            self.get_lc_ranges()
            self.process_lowcomplex( variants_proc )
            self.variants = variants_proc

        # process indel removal filter
        if self.indel:
            variants_proc = {}
            self.process_indel( variants_proc )
            self.variants = variants_proc

        # run filters on data in INFO field
        if self.mq or self.dp:
            variants_proc = {}
            self.process_info_field( variants_proc )
            self.variants = variants_proc

        # run filters on data in genotype field
        if self.gq or self.gt:
            variants_proc = {}
            self.process_genotype_field( variants_proc )
            self.variants = variants_proc


################
##### MAIN #####
################
def main( prog_name, argv ):
    # ARG PROCESSING
    parser = argparse.ArgumentParser( prog=prog_name, description='process out SNPs from vcf located in ranges specified by gff file, only filters specified will be processed',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument('-i','--vcf_file', dest='vcf_file', metavar='VCF_FILE', type=str, required=True, help='file containing SNPs and other variants to be processed')
    parser.add_argument('-g','--gff_file', dest='gff_file', metavar='GFF_FILE', type=str, help='annotation file containing the list of ranges specified as low complexity by RepeatMasker')
    parser.add_argument('-z','--indel', dest='indel', action='store_true', default=False, help='indel removal switch')
    parser.add_argument('-m','--mq', dest='mq', metavar='MQ', type=float, help='minimum mapping quality filter')
    parser.add_argument('-d','--dp', dest='dp', metavar='DP', type=float, help='minimum DP value filter')
    parser.add_argument('-q','--gq', dest='gq', metavar='GQ', type=float, help='minimum GQ value filter')
    parser.add_argument('-t','--gt', dest='gt', metavar='GT', type=str, action='append', help='genotype to filter out in gt format( eg, 0/1 )')
    parser.add_argument('-n','--no_count', dest='no_count', action='store_true', default=False, help='suppress removal counting and its output')
    parser.add_argument('-o','--outfile', dest='outfile', metavar='OUTFILE', type=str, help='output file to write to')
    args = parser.parse_args(argv)

    # SET VARIBLES FROM COMMAND LINE ARGS
    vcf_file = args.vcf_file
    gff_file = args.gff_file
    indel = args.indel
    mq = args.mq
    dp = args.dp
    gq = args.gq
    gt = args.gt
    outfile = args.outfile
    no_count = args.no_count
  
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

    # initialize Filter object
    fltr = Filter( gff_file, indel, mq, dp, gq, gt, no_count )

    # READ VARIANTS
    fltr.get_variants( vcf_file )

    # PROCESS FILTERS
    fltr.run_filters()

    # WRITE FILTERED VARIANTS
    fltr.write_snps( outfile )

    # WRITE COUNT REMOVALS
    if not no_count:
        fltr.write_rem_counts()


if __name__=='__main__':
    main(sys.argv[0], sys.argv[1:])
