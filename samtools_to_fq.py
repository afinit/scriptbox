#!/usr/bin/env python
# create fastq file from output of samtools

import sys

filename = sys.argv[1]
fh = open( filename, 'w' )

for line in sys.stdin:
  line = line.split()
  fh.write( '@' + line[0] + '\n' )
  fh.write( line[9] + '\n' )
  fh.write( '+\n' )
  fh.write( line[10] + '\n' )

fh.close()

