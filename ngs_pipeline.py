from ruffus import *
import argparse
from logbook import Logger, FileHandler
import getpass
import subprocess
import sys
import datetime



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

i = datetime.datetime.now()
date_time_iso = '%s' %i.isoformat()

print(date_time_iso)