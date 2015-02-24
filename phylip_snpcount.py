#!/usr/bin/env python
# reads in phylip formatted files where each genome is unbroken on one line, then it does pairwise comparisons of genomes to come up with total SNPs between the each genome pair
#   output is a 2d matrix in a file
# USAGE: phylip_snpcount.py -p file.phy [-o file.txt]

import os, sys
import argparse


# read in variants from vcf file
def read_phylip_file(phylip_file, genomes, genome_names):
    print 'READING {0}'.format( phylip_file )
    with open( phylip_file, 'r' ) as fh:
        line = fh.readline().split()
        try:
            int( line[0] )
            int( line[1] )
        except ValueError:
            print 'Error: improper format: PHYLIP files should have two integers in the first line describing how many genomes are present and the length of the sequences'
            sys.exit()

        for line in fh:
            line = line.split()
            
            if len(line) != 2:
                print 'Error: Line does not have name and sequence'
                sys.exit()

            # add genome to genomes list
            genome_names.append(line[0])
            genomes.append(line[1])


################
##### MAIN #####
################
def main( prog_name, argv ):
    # ARG PROCESSING
    parser = argparse.ArgumentParser( prog=prog_name, description='reads in and stores multiple vcf files',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument('-p','--phylip', dest='phylip_file', metavar='PHYLIP_FILE', type=str, required=True, help='phylip formatted file to count SNPs in')
    parser.add_argument('-o','--outfile', dest='outfile', metavar='OUTFILE', type=str, help='file to print matrix to')
    args = parser.parse_args(argv)

    # SET VARIBLES FROM COMMAND LINE ARGS
    phylip_file = args.phylip_file
    outfile = args.outfile

    if outfile is None:
        outfile = os.path.splitext( phylip_file )[0] + '_snpcount.txt'
  
    # check if input files exist
    if not os.path.isfile( phylip_file ):
        print 'Error: {0} does not exist'.format( ref_file )
        sys.exit(1)

    # INITIALIZE
    genomes = []
    genome_names = []
 
    # get genomes from phylip file
    read_phylip_file( phylip_file, genomes, genome_names )
    
    # initialize snp_count matrix
    snp_count = [[0 for i in genomes] for j in genomes]
    total_genomes = len(genomes)
    seq_len = len(genomes[0])

    # process snps
    print "PROCESSING SNPS"
    for pos in xrange(seq_len):
        base_same = set([x for x in xrange(total_genomes)])
        
        # loop until no more differences or matches are found
        while(1):
            first_index = list(base_same)[0]

            # find indices different from the first in the list
            base_diff = set([x for x in base_same if genomes[first_index][pos] != genomes[x][pos]])
            
            # if base_diff is empty, break the loop
            if base_diff == set():
                break

            # remove indices of different bases from base_same
            base_same = base_same - base_diff

            # loop through base_same and base_diff to increase all of the appropriate counters
            for i in base_same:
                for j in base_diff:
                    snp_count[i][j] += 1
                    snp_count[j][i] += 1
            
            # set base_same to the base_diff set to process those indices
            base_same = base_diff



    # print matrix
    for i in xrange(total_genomes):
        print genome_names[i] + '\t' + str(max(snp_count[i])) + '\t' + str(sum(snp_count[i])/float(seq_len)) 

    print 'WRITING SNPCOUNT MATRIX TO FILE'
    with open( outfile, 'w' ) as out_fh:
        out_fh.write( '#gen\t' + '\t'.join( genome_names ) + '\n' )
        for i in xrange( total_genomes ):
            out_fh.write( genome_names[i] + '\t' + '\t'.join( str(n) for n in snp_count[i] ) + '\n' )

    print 'Done'


if __name__=='__main__':
    main(sys.argv[0], sys.argv[1:])
