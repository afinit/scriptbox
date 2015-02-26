#!/usr/bin/env python
# input files:  read_map.bam
#               ref.fa
#   ref.fa must be indexed using samtools faidx and picardtools CreateSequenceDictionary.jar
# usage: snp_call.py read_map.bam ref.fa

import sys, os
import subprocess

def print_usage():
    print "Usage: snp_call.py read_map.bam ref.fa [start_cmd]"
    print ' 1. MarkDuplicates'
    print ' 2. RealignerTargetCreator'
    print ' 3. IndelRealigner'
    print ' 4. HaplotypeCaller'
    print
    print '   Note: ref.fa must be indexed using samtools faidx and picardtools CreateSequenceDictionary.jar'
    sys.exit()

def build_bam_index( map_file ):
  try:
    print "Starting bam index"
    subprocess.check_call( 'java -Xmx4g -jar ~/src/picard-tools*/BuildBamIndex.jar INPUT= ' + map_file, shell=True )
  except subprocess.CalledProcessError as e:
    print 'Error: failed running picard-tools BuildBamIndex: ' + str(e.returncode)
    print_usage()


def main( argv ):
  start_cmd = 0

  # set tmp dir in user home directory
  tmp_dir = os.path.expanduser('~') + os.sep + 'tmp'
  if not os.path.isdir( tmp_dir ):
    print 'Creating tmp directory: ' + tmp_dir
    os.mkdir( tmp_dir )

  if len( argv ) < 2:
    print "Error: arguments missing"
    print_usage()
  elif len( argv ) > 2:
    start_cmd = int(argv[2])

  map_file = argv[0]
  map_file_base = os.path.splitext(map_file)[0]
  ref_file = argv[1]
  ref_file_base = os.path.splitext(ref_file)[0]
  exist_flag = True
  output_file = map_file_base

  # check for existence of mapping file and reference file
  if not os.path.isfile(map_file) or not os.path.splitext(map_file)[1] == '.bam':
    print "Error: Can't find read mapping file: " + map_file
    exist_flag = False
  if not os.path.isfile(ref_file) or not os.path.splitext(ref_file)[1] == '.fa':
    print "Error: Can't find reference file: " + ref_file
    exist_flag = False

  if not exist_flag:
    print_usage()
  
  # check for dict index
  if os.path.isfile(ref_file_base + '.dict'):
    if os.path.getctime(ref_file_base + '.dict') < os.path.getctime(ref_file):
      os.remove(ref_file_base + '.dict')
    
  if not os.path.isfile(ref_file_base + '.dict'):
    print "Building dict index: " + ref_file_base + '.dict'
    try:
      subprocess.check_call( 'java -Djava.io.tmpdir=' + tmp_dir + ' -jar ~/src/picard-tools*/CreateSequenceDictionary.jar R= ' + ref_file + ' O= ' + ref_file_base + '.dict', shell=True )
    except subprocess.CalledProcessError as e:
      print 'Error: failed running picard-tools CreateSequenceDictionary: ' + str(e.returncode)
      print_usage()

  # check for faidx index
  if os.path.isfile(ref_file + '.fai'):
    if os.path.getctime(ref_file + '.fai') < os.path.getctime(ref_file):
      os.remove(ref_file + '.fai')

  if not os.path.isfile(ref_file + '.fai'):
    print "Building faidx index: " + ref_file + '.fai'
    try:
      subprocess.check_call( 'samtools faidx ' + ref_file, shell=True )
    except subprocess.CalledProcessError as e:
      print 'Error: failed running samtools faidx: ' + str(e.returncode)
      print_usage()

  if start_cmd <= 1:
    # build bam index picard tools
    build_bam_index( map_file )

    # MarkDuplicates
    try:
      print 'MarkDuplicates'
      subprocess.check_call( 'java -Djava.io.tmpdir=' + tmp_dir + ' -Xmx4g -jar ~/src/picard-tools*/MarkDuplicates.jar INPUT=' + map_file_base + '.bam OUTPUT=' + map_file_base + '.dedup.bam METRICS_FILE=' + map_file_base + '.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000', shell=True )
    except subprocess.CalledProcessError as e:
      print 'Error: failed running picard-tools MarkDuplicates: ' + str(e.returncode)
      print_usage()

  if start_cmd <= 2:
    # build bam index picard tools: dedup
    build_bam_index( map_file_base + '.dedup.bam' )

    # RealignerTargetCreator
    try:
      print 'RealignerTargetCreator'
      subprocess.check_call( 'java -Djava.io.tmpdir=' + tmp_dir + ' -Xmx4g -jar ~/src/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ' + ref_file + ' -o ' + map_file_base + '.list -I ' + map_file_base + '.dedup.bam', shell=True )
    except subprocess.CalledProcessError as e:
      print 'Error: failed running GATK RealignerTargetCreator: ' + str(e.returncode)
      print_usage()
   
  if start_cmd <= 3:
    # IndelRealigner
    try:
      print 'IndelRealigner'
      subprocess.check_call( 'java -Djava.io.tmpdir=' + tmp_dir + ' -Xmx4g -jar ~/src/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ' + ref_file + ' -targetIntervals ' + map_file_base + '.list -I ' + map_file_base + '.dedup.bam -o ' + map_file_base + '.dedup.realign.bam', shell=True )
    except subprocess.CalledProcessError as e:
      print 'Error: failed running GATK IndelRealigner: ' + str(e.returncode)
      print_usage()
  
  # build bam index picard tools: realign
  build_bam_index( map_file_base + '.dedup.realign.bam' )

  # HaplotypeCaller
  try:
    print 'HaplotypeCaller'
    subprocess.check_call( 'java -Djava.io.tmpdir=' + tmp_dir + ' -Xmx4g -jar ~/src/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ' + ref_file + ' -l INFO -I ' + map_file_base + '.dedup.realign.bam -o ' + output_file + '.vcf', shell=True )
  except subprocess.CalledProcessError as e:
    print 'Error: failed running GATK HaplotypeCaller: ' + str(e.returncode)
    print_usage()
  
if __name__=='__main__':
  main(sys.argv[1:])
