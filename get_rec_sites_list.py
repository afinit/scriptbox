#!/usr/bin/env python
# pulls data from neb site, parses and saves a tabular file of neb recognition sites that do not contain ambiguous nucleotides

import urllib2
import re

wbpg = urllib2.urlopen( 'https://www.neb.com/tools-and-resources/selection-charts/alphabetized-list-of-recognition-specificities' ).read()

wbpg = wbpg.split( '<table>' )[1]
wbpg = wbpg.split( '</table>' )[0]

# remove the table head section
wbpg = wbpg.split( '<tr>', 2 )[2]

wbpg_l = wbpg.split( '</tr><tr>' )

site_list = dict()

for site in wbpg_l:
  seq = ''
  seq_l = ''
  site_id = ''
  rem = []

  site_l = re.findall( '<td>(.*?)</td>', site )
  #process the sequence first
  seq_l = list( site_l[0] )
  for c in seq_l:
    if c.isalpha() and c not in ['A','T','C','G']:
      seq_l = ''
      break
    elif c not in ['A','T','C','G']:
      #build list of characters to remove from string
      rem.append( c )

  # bail if ambiguous nucleotides are found
  if seq_l == '':
    continue

  # remove non nucleotide sequences from seq
  for c in rem:
    seq_l.remove( c )

  seq = "".join( seq_l )
  
  #process the recognition site name second
  site_id = re.findall( '<a.*?>(.*?)</a>', site_l[1] )[0]
  site_id = re.sub( '<.*>', '', site_id )
  
  site_list[ site_id ] = seq

fh = open( 'rec_sites.dat', 'w' )

for key in site_list.keys():
  fh.write( key + '\t' + site_list[key] + '\n' )

fh.close()
