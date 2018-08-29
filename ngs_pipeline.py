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

########################################################

# Second quality check- FASTQC

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


########################################################

# Aligning to a reference genome- BWA

# @follows(second_quality_check)
# @collate('{input}*.qfilter.fastq.gz'.format(input=args.input_dir), formatter('([^/]+)_L001_R[12]_001.qfilter.fastq.gz$'),'{path[0]}/{1[0]}.bwa.bam')
# def align_FASTQ(infile,outfile):

# 	FASTQ1 = infile[0]
# 	FASTQ2 = infile[1]

# 	# Removes the qfilter from fastq file name
# 	original_FASTQ = re.sub(r'.qfilter.fastq.gz','.fastq.gz',FASTQ1)
# 	# Ensures renamed file is in the correct location
# 	path = os.path.realpath('%s' % original_FASTQ)
# 	sample_name = re.sub(r'_L001_R1_001.qfilter.fastq.gz','',FASTQ1)

# 	try:
# 		# Extracting the patient ID for the storage of the bam file
# 		flowcell = re.search(r'/([0-9]{6})_(M.+?)_([0-9].+?)-(.+?)/',path).groups(1)[3]
# 		split_path = path.split('/')

# 	except AttributeError:

# 		flowcell = 'NA'

# 	align_command = ('{bwa} mem -t {threads} -M -k 18 -R "@RG\\tID:{worksheet}.{flowcell}\\tCN:WMRGL\\tDS:{panel}\\tDT:{date}\\tSM:{sample_name}\\tLB:{worksheet}\\tPL:ILLUMINA" \
# 				{genome} {FASTQ1} {FASTQ2} \
# 	 			| sed \'s/-R @RG.*//\' - \
# 	 			| {samtools} view -Sb - \
# 	 			| {samtools} sort -T {sample_name}.temp -O bam - > {outfile}'.format(

# 				samtools    = config_dict['samtools'], 
# 				bwa         = config_dict['bwa'],
# 				threads     = config_dict['bwa_threads'],
# 				worksheet   = 'test',
# 				flowcell    = flowcell,
# 				date        = date_time_iso,
# 				sample_name = sample_name,
# 				genome      = config_dict['reference_genome'],
# 				FASTQ1      = FASTQ1,
# 				FASTQ2      = FASTQ2,
# 				outfile     = outfile, 
# 				panel       ='panel'))

# 	#Log information
# 	log.info('aligning: '+ FASTQ1 + ' ' + FASTQ2)
# 	log.info(align_command)
# 	#Run command
# 	os.system(align_command)


#################################################################

# Index the bam file (Sorting completed in the previous step)
# Indexing allows fast random accessing

# @follows(align_FASTQ)
# @transform('{input}*.bwa.bam'.format(input=args.input_dir),suffix('.bwa.bam'),'.bwa.bam.bai')
# def index_original_bam(infile,outfile):
# 	command = '{samtools} index {infile}'.format(
# 	samtools = config_dict['samtools'],
# 	infile = infile
# 	)

# 	log.info('indexing file : '+ infile)
# 	log.info(command)
# 	os.system(command)


#################################################################
# Local re-alignment- Abra (version 1 (output needs sorting))

# @follows(index_original_bam)
# @transform('{input}*.bwa.bam'.format(input=args.input_dir),suffix('.bwa.bam'),'.bwa.realn.bam')
# def abra_realign_bam(infile,outfile):
# 	# make a temporary directory for ABRA to work in:
# 	temp_folder = '{dir}/temp'.format(dir=args.input_dir)

# 	if os.path.isdir(temp_folder) != True:
# 		os.mkdir(temp_folder)


# 	name = re.sub(".bwa.bam", "", infile)
# 	command = "java -Xmx4G -jar {abra}\
# 	--in {infile} \
# 	--out {outfile} \
# 	--ref {ref} \
# 	--targets {targets}\
# 	--threads {threads} \
# 	--working {working} > {log_location} 2>&1".format(

# 		abra=config_dict['abra'],
# 		infile=infile,
# 		outfile=outfile,
# 		ref=config_dict['reference_genome'],
# 		targets=config_dict['bed_file'],
# 		threads=config_dict['abra_threads'],
# 		working=temp_folder,
# 		log_location=name+'.abra.log'
# 		)

# 	os.system(command)



#################################################################

# Sorting re-aligned bam file

# @follows(abra_realign_bam)
# @transform(["{input}*.bwa.realn.bam".format(input=args.input_dir)], suffix(".bwa.realn.bam"), ".bwa.realn.sorted.bam")
# def sorted_realign_bam(infile, outfile):



# 	#logging
# 	log.info('sorting realigned bam file: ' + infile)
# 	command = '{samtools} sort {infile} > {outfile}'.format(
# 		samtools=config_dict['samtools'],
# 		infile=infile,
# 		outfile=outfile)
# 	log.info(command)
# 	#run command
# 	os.system(command)


################################################################

# Indexing sorted-realigned bam file


# @follows(sorted_realign_bam)
# @transform(['{input}*.bwa.realn.sorted.bam'.format(input=args.input_dir)], suffix(".bwa.realn.sorted.bam"), ".bwa.realn.sorted.bam.bai" )
# def index_realigned_bam(infile,outfile):
# 	command = '{samtools} index {infile}'.format(
# 	samtools = config_dict['samtools'],
# 	infile = infile
# 	)

# 	log.info('indexing file : '+ infile)
# 	log.info(command)
# 	os.system(command)

################################################################

# Fixing mates

# @follows(index_realigned_bam)
# @transform(['{input}*.bwa.realn.sorted.bam'.format(input=args.input_dir)], suffix('bwa.realn.sorted.bam'), 'bwa.realn.fixed.bam')
# def picard_fix_mate(infile,outfile):

# 	picard_temp_folder = '{dir}/picar_temp_folder'.format(dir=args.input_dir)

# 	if os.path.isdir(picard_temp_folder) != True:
# 		os.mkdir(picard_temp_folder)

# 	command = 'java -Xmx2G -jar {picard} FixMateInformation INPUT={infile} OUTPUT={outfile} TMP_DIR={picard_temp_dir}'.format(
# 		picard = config_dict['picard'],
# 		infile = infile,
# 		outfile = outfile,
# 		picard_temp_dir = picard_temp_folder
# 		) 

# 	log.info('Mate fixing: '+ infile)
#  	log.info(command)
# 	os.system(command)


################################################################

# Indexing fixed bam file


# @follows(picard_fix_mate)
# @transform(['{input}*.bwa.realn.fixed.bam'.format(input=args.input_dir)], suffix(".bwa.realn.fixed.bam"), ".bwa.realn.fixed.bam.bai" )
# def index_fixed_bam(infile,outfile):
# 	command = '{samtools} index {infile}'.format(
# 	samtools = config_dict['samtools'],
# 	infile = infile
# 	)

# 	log.info('indexing fixed file: '+ infile)
# 	log.info(command)
# 	os.system(command)

####################################################################

# Sorting fixed bam file


# @follows(index_fixed_bam)
# @transform(["{input}*.bwa.realn.fixed.bam".format(input=args.input_dir)], suffix(".bwa.realn.fixed.bam"), ".bwa.realn.fixed.sorted.bam")
# def sorted_fixed_bam(infile, outfile):
# 	command = '{samtools} sort {infile} > {outfile}'.format(
# 		samtools=config_dict['samtools'],
# 		infile=infile,
# 		outfile=outfile)
# 	log.info(command)
# 	#run command
# 	os.system(command)

#####################################################################

# Indexing the realigned, fixed and sorted bam
# @follows(sorted_fixed_bam)
# @transform(["{input}*.bwa.realn.fixed.sorted.bam".format(input=args.input_dir)], suffix('.bwa.realn.fixed.sorted.bam'), '.bwa.realn.fixed.sorted.bam.bai') 
# def index_fixed_sorted_bam(infile,outfile):
# 	command = '{samtools} index {infile}'.format(
# 	samtools = config_dict['samtools'],
# 	infile = infile
# 	)

# 	log.info('indexing fixed sorted file: '+ infile)
# 	log.info(command)
# 	os.system(command)


#####################################################################

# Variant calling- VarScan

####### VarScan must be run in germline mode because there are no control samples to run #######

# @follows(index_fixed_sorted_bam)
# @transform(["{input}*.bwa.realn.fixed.sorted.bam".format(input=args.input_dir)], suffix('.bwa.realn.fixed.sorted.bam'), '.vcf') 
# def call_variants_varscan(infile,outfile):

# 	name = re.sub(r'.bwa.realn.fixed.sorted.bam', '', infile)
# 	samp_name = name.split('/')[-1]

	## Call all snp variants and store as VCF
	# snp_command = "{samtools} mpileup -B \
	# 		-f {reference} \
	# 		{infile} | java -jar {var_scan} mpileup2snp \
	# 		--min-coverage 8 \
	# 		--min-reads2 2 \
	# 		--min-avg-qual 10 \
	# 		--p-value 99e02 \
	# 		--output-vcf 1 \
	# 		--strand-filter 0 \
	# 		> {name}.snps.vs2.vcf".format(samtools=config_dict['samtools'],reference=config_dict['reference_genome'], infile=infile, var_scan=config_dict['varscan'], name=name)
	# log.info('Calling variants for {file}'.format(file=infile))
	# log.info(snp_command)
	# os.system(snp_command)


	## Call all indel variants and store as VCF
	# indel_command = "{samtools} mpileup -B \
	# 		-f {reference} \
	# 		{infile} | java -jar {var_scan} mpileup2indel \
	# 		--min-coverage 8 \
	# 		--min-reads2 2 \
	# 		--min-avg-qual 10 \
	# 		--min-avg-freq .01 \
	# 		--p-value 99e-02 \
	# 		--output-vcf 1 \
	# 		--strand-filter 0 \
	# 		> {name}.indel.vs2.vcf".format(samtools=config_dict['samtools'], reference=config_dict['reference_genome'], infile=infile, var_scan=config_dict['varscan'], name=name)

	# log.info('Calling indel variants for {file}'.format(file=samp_name))
	# log.info(indel_command)
	# os.system(indel_command)


	# # Must add perl5lib to path before vcftools can be used
	# # replace the 'Sample1' label
	# sed_snp_command = "sed -i 's/Sample1/{samp_name}/' {name}.snps.vs2.vcf".format(samp_name=samp_name, name=name)
	# sed_indel_command = "sed -i 's/Sample1/{samp_name}/' {name}.indel.vs2.vcf".format(samp_name=samp_name, name=name)
	# os.system(sed_snp_command)
	# os.system(sed_indel_command)
	# log.info('Removing Sample1 from the vcf files')

	# combine_vcf_command = "{vcf_concat} {name}.snps.vs2.vcf {name}.indel.vs2.vcf > {name}.merged.vcf".format(vcf_concat=config_dict['vcf_concat'], name=name)
	# os.system(combine_vcf_command)
	# log.info('Combining SNP and indel VCF files')

##################################################################

## variant normalization- VT
## Splits, left-align and trims in one step

# @follows(call_variants_varscan)
# @transform(['{input}*.merged.vcf'.format(input=args.input_dir)], suffix('.merged.vcf'), '.normal.vcf')
# def normalize_vcf(infile,outfile):
# 	command = '{vt} decompose -s {sample} | {vt} normalize -m -r {reference} -o {output} -'.format(
# 		vt=config_dict['vt'],
# 		sample=infile,
# 		reference=config_dict['reference_genome'],
# 		output=outfile)

# 	os.system(command)
# 	log.info('Decomposing and norimalizing merged VCF file')
# 	log.info(command)


# @follows(normalize_vcf)
@transform(['{input}*.normal.vcf'.format(input=args.input_dir)], suffix('.normal.vcf'), '.vep.vcf')
def vep_annotation(infile,outfile):
	command = '{vep} -i {infile} -o {outfile} --cache --fork 4 --dir_cache=/media/sf_S_drive/ensembl/cache/ --vcf --flag_pick \
			--exclude_predicted --everything --dont_skip --total_length --offline --fasta {reference}'.format(
				vep= config_dict['vep'],
				infile= infile,
				outfile= outfile,
				reference= config_dict['reference_genome']
				)

	os.system(command)

pipeline_run()