#!/usr/bin/env python

import argparse
import sys



#-------------------------------------------------------------------



def checkDirectories(params):
	error = False
	dirs = [params['SAMPLES_DIR'], params['CONFIG_DIR'], params['RESULTS']\
		, params['TAXONOMY_FILE'], params['SCRIPTS'], params['BINARIES'], params['PUFFERFISH_INDEX']]

	if '' in dirs:
		error = True
	if '""' in dirs:
		error = True
	if params['MODE'] == '"2"':
		if params['HOST_INDEX'] == '' or params['HOST_INDEX'] == '""':
			error = True

	return error



def checkParameters(params):
	error = False
	if params['MODE'] not in ('"1"', '"2"'):
		 error = True
	if params['STRATEGY'].lower() not in ('"pe"', '"se"'):
		error = True
	if params['ALGORITHM'].lower() not in ('"agamemnon"', '"kallisto"', '"kraken"'):
		error = True
	if params['HOST_SAM'].lower() not in ('"true"', '"false"'):
		error = True
	if params['TO_BAM'].lower() not in ('"true"', '"false"'):
		error = True
	if params['CLEAR_ALL'].lower() not in ('"true"', '"false"'):
		error = True
	if params['LIBTYPE'].lower() == '""':
		error = True

	return error



def readParameters(configF):
	params = dict()
	with open(configF) as f:
		for line in f:
			if line.startswith('#') or not line.strip():
				continue
			else:
				elements = line.split(':')
				params[elements[0].strip()] = elements[1].strip()

	return params



def writeOutput(error, configValidation):
	configValidation.write('%s\n' % (str(error)))



#-------------------------------------------------------------------



if __name__ == '__main__':
	parser = argparse.ArgumentParser(__file__, description = 'Configuration file check')
	parser.add_argument('-configDir', '-cd', help = 'Config file directory', default = None)
	parser.add_argument('-outDir', '-od', help = 'Output directory', default = None)

	args = parser.parse_args()
	configValidation = open(args.outDir + '/configuration_file.report', 'w')
	params = readParameters(args.configDir)
	error = checkDirectories(params)
	if error:
		writeOutput(error, configValidation)
	else:
		error = checkParameters(params)
		writeOutput(error, configValidation)




