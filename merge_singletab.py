#!/usr/bin/env python
# this will accept an optional key_file, key_field, value_field, and outfile, it merges multiple files that contain lists of tab_delimited files in similar formats using the key_field as the first field of the new file and the value_field is used as the value 
# USAGE: merge_singletab.py [-f <key_file.txt>] [-k <key_field>] [-v <value_field] [-o <outfile>] <FILES_TO_BE_MERGED>

import os, sys
import argparse


def process_file( infile, key_field, value_field, key_list=set() ):
  data_dict = {}
  
  if len(key_list) == 0:
    add_keys = True
  else:
    # initialize data dictionary
    for i in key_list:
      data_dict[i] = ''
    add_keys = False

  with open( infile, 'r' ) as fh:
    # get first line
    line = fh.readline().split()

    # loop through file 
    while line != []:
      # prevent index error
      if len(line) <= key_field or len(line) <= value_field:
        print 'Error: Not enough fields in {}'.format( infile )
        sys.exit()

      # if key_list is provided to this function, check to see if the current line's key_field is in key_list
      if add_keys:
        data_dict[line[key_field]] = line[value_field]
      else:
        if line[key_field] in key_list:
          data_dict[line[key_field]] = line[value_field]

      line = fh.readline().split()

  return data_dict

def main( prog_name, argv ):
  # ARG PROCESSING
  parser = argparse.ArgumentParser( prog=prog_name, description='merges multiple files that contain lists of tab_delimited files in similar formats using the key_field as the first field of the new file and the value_field is used as the value',
      formatter_class=argparse.ArgumentDefaultsHelpFormatter )
  parser.add_argument('infiles', metavar='INFILES', type=str, nargs='+', help='files to be merged')
  parser.add_argument('-f,--key_file', dest='key_file', metavar='KEY_FILE', type=str, help='file containing the list of keys to use in the processing')
  parser.add_argument('-k,--key_field', dest='key_field', metavar='KEY_FIELD', type=int, default=0, help='field in the INFILES to use as key')
  parser.add_argument('-v,--value_field', dest='value_field', metavar='VALUE_FIELD', type=int, default=1, help='field in the INFILES to use as value')
  parser.add_argument('-o,--outfile', dest='outfile', metavar='OUTFILE', type=str, default='merged_out.txt', help='output file to which the merged data will be written')
  args = parser.parse_args(argv)

  # SET VARIBLES FROM COMMAND LINE ARGS
  infiles = args.infiles
  files_to_process = [x for x in infiles]
  key_file = args.key_file
  key_field = args.key_field
  value_field = args.value_field
  outfile = args.outfile
  data_values = {}

  # if a key_file is specified, read this file here
  if not key_file is None:
    with open( key_file, 'r' ) as key_fh:
      key_list = set()

      # initialize line holding variable
      line = key_fh.readline().rstrip()

      # loop through gene_file
      while line != '':
        key_list.add( line )
        line = key_fh.readline().rstrip()
    
    # build dictionary for data values
    for i in key_list:
      data_values[i] = []

  else:
    data_dict = process_file( files_to_process.pop(0), key_field, value_field )
    key_list = data_dict.keys()
    
    # build dictionary for data values
    for i in key_list:
      data_values[i] = []
      data_values[i].append( data_dict[i] )
    
  # loop through file list to pull out data
  for f in files_to_process:
    data_dict = process_file( f, key_field, value_field, key_list )
    
    for i in key_list:
      data_values[i].append( data_dict[i] )

  # print to csv file
  with open( outfile, 'w' ) as out_fh:
    out_fh.write( 'keys' )
    for i in infiles:
      out_fh.write( ', ' + os.path.splitext(os.path.basename(i))[0] )

    out_fh.write( '\n' )

    for i in data_values:
      out_fh.write( i + ', ' + ', '.join( data_values[i] ) + '\n' )


if __name__=='__main__':
  main(sys.argv[0], sys.argv[1:])
