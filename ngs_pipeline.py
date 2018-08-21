from ruffus import *
import argparse
from logbook import Logger, FileHandler
import getpass
import subprocess
import sys
import datetime
import glob
import os
import re


###################################################
# Open config file and extract paths to tools
# Add to a config dictionary

config_dict = {}
with open('config.txt', 'r') as config:
	for line in config:
		line = line.strip().split('=', 1)
		config_dict[line[0]] = line[1]

###################################################
# Adding parse arguements

parser = argparse.ArgumentParser(description='An NGS pipeline')
parser.add_argument('--input_dir', nargs='?')
parser.add_argument('--config', nargs='?')
args = parser.parse_args()


####################################################
# Create a logger

FileHandler('{input_dir}ngs_pipeline.log'.format(input_dir=args.input_dir)).push_application()
log = Logger('pipeline_logger')
log.info('Starting the pipeline as user ' + getpass.getuser())
git_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD']).strip()
log.info('pipeline git hash: '+ git_hash)
log.info('input_dir: '+ args.input_dir)
log.info('config: '+ args.config)
log.info('using python version: '+ sys.version.strip())
log.info('fastqc version: '+ subprocess.check_output([config_dict['fastqc'], '--version']).strip())

####################################################

# Set current time

i = datetime.datetime.now()
date_time_iso = '%s' %i.isoformat()

####################################################

# Get a list of input files

# Search for all fastq.gz files in the input directory
initial_files = glob.glob('{input_dir}*.fastq.gz'.format(input_dir=args.input_dir))

####################################################

# Checking quality- FASTQC
# Create QC reports folder and conduct quality analysis

# @mkdir('{input}qc_reports'.format(input=args.input_dir))
# @transform(initial_files, suffix('.fastq.gz'), ' ')
# def first_quality_check(infile, outfile):
# 	fastqc_cmd = '{fastqc} {filename} --outdir={output_dir}'.format(
# 		fastqc=config_dict['fastqc'],
# 		filename=infile,
# 		output_dir='{input}qc_reports'.format(input=args.input_dir)
# 		)
# 	log.info('first qc check for file: '+ infile)
# 	log.info(fastqc_cmd)
# 	os.system(fastqc_cmd)

####################################################

# Trimming low quality reads- Trimmomatic

# @follows(first_quality_check)
# path[0] the whole path to the whole directory
# 1[0] refers to the group 1 (first group) in the formatter regex (group 0 would be the whole match) 
# @collate(initial_files, formatter('([^/]+)R[12]_001.fastq.gz$'), '{path[0]}/{1[0]}.fastq.gz')
# def trim_low_quality_reads(infile, outfile):

# 	'''
# 	Trimmomatic command:
# 	See documentaiton for parameters used
# 	2 input files- FASTQ1 FASTQ2
# 	4 output files- 2 (paired) qfilter files
# 					2 unpaired files (where a read survived, but a partener did not)
# 	'''

# 	FASTQ1 = infile[0]
# 	FASTQ2 = infile[1]

#  	trim_cmd = 'java -jar {trimmomatic} PE -phred33 {R1} {R2} {filename1}.qfilter.fastq.gz \
#  	{filename1}.unpaired.fastq.gz {filename2}.qfilter.fastq.gz {filename2}.unpaired.fastq.gz \
#  	ILLUMINACLIP:{illuminaclip} CROP:150 SLIDINGWINDOW:4:20 MINLEN:50'.format(
#  	trimmomatic=config_dict['trimmomatic'],
#  	R1=FASTQ1,
#  	R2=FASTQ2,
#  	filename1=FASTQ1[:-9],
#  	filename2=FASTQ2[:-9],
#  	illuminaclip=config_dict['illuminaclip']
#  	)

#  	log.info('quality trim for files: '+ FASTQ1 + ' ' + FASTQ2)
#  	log.info(trim_cmd)
#  	os.system(trim_cmd)

# @follows(quality_trim)
# @transform('{input}*.qfilter.fastq.gz'.format(input=args.input_dir), suffix('.qfilter.fastq.gz'), ' ')
# def second_quality_check(infile, outfile):
# 	fastqc_cmd = '{fastqc} {filename} --outdir={output_dir}'.format(
# 		fastqc=config_dict['fastqc'],
# 		filename=infile,
# 		output_dir='{input}qc_reports'.format(input=args.input_dir)
# 		)

# 	log.info('second qc check for file: '+ infile)
# 	log.info(fastqc_cmd)
# 	os.system(fastqc_cmd)


# @follows(second_quality_check)
@collate('{input}*.qfilter.fastq.gz'.format(input=args.input_dir), formatter('([^/]+)_L001_R[12]_001.qfilter.fastq.gz$'),'{path[0]}/{1[0]}.bwa.bam')
def align_FASTQ(infile,outfile):

	FASTQ1 = infile[0]
	FASTQ2 = infile[1]

	# Removes the qfilter from fastq file name
	original_FASTQ = re.sub(r'.qfilter.fastq.gz','.fastq.gz',FASTQ1)
	# Ensures renamed file is in the correct location
	path = os.path.realpath('%s' % original_FASTQ)
	sample_name = re.sub(r'_L001_R1_001.qfilter.fastq.gz','',FASTQ1)

	try:
		# Extracting the patient ID for the storage of the bam file
		flowcell = re.search(r'/([0-9]{6})_(M.+?)_([0-9].+?)-(.+?)/',path).groups(1)[3]
		split_path = path.split('/')

	except AttributeError:

		flowcell = 'NA'

	align_command = ('{bwa} mem -t {threads} -M -k 18 -R "@RG\\tID:{worksheet}.{flowcell}\\tCN:WMRGL\\tDS:{panel}\\tDT:{date}\\tSM:{sample_name}\\tLB:{worksheet}\\tPL:ILLUMINA" \
				{genome} {FASTQ1} {FASTQ2} \
	 			| sed \'s/-R @RG.*//\' - \
	 			| {samtools} view -Sb - \
	 			| {samtools} sort -T {sample_name}.temp -O bam - > {outfile}'.format(

				samtools    = config_dict['samtools'], 
				bwa         = config_dict['bwa'],
				threads     = config_dict['bwa_threads'],
				worksheet   = 'test',
				flowcell    = flowcell,
				date        = date_time_iso,
				sample_name = sample_name,
				genome      = config_dict['reference_genome'],
				FASTQ1      = FASTQ1,
				FASTQ2      = FASTQ2,
				outfile     = outfile, 
				panel       ='panel'))

	#Log information
	log.info('aligning: '+ FASTQ1 + ' ' + FASTQ2)

	log.info(align_command)

	#Run command
	os.system(align_command)

pipeline_run()