#!/usr/bin/env python3

import argparse
import os
import yaml



def read_size(sizeDir):
	size = list()
	with open(sizeDir + '/indexSize.tmp') as f:
		for line in f:
			size.append(line.strip())

	return size[0]



def edit_yaml(yamlDir, size):
	yamlFile = yamlDir + '/config.yml'
	with open(yamlFile) as f:
		values = yaml.safe_load(f)
	values['resources']['MEM_MB'] = int(size)
	with open(yamlFile, 'w') as f:
		yaml.safe_dump(values, f, default_flow_style = False)

	return None




if __name__ == '__main__':
	parser = argparse.ArgumentParser(__file__, description = "Edit config yaml file")
	parser.add_argument('-yamlDir', '-yd', help = 'Yaml file directory', default = None)
	parser.add_argument('-sizeDir', '-sd', help = 'Index size file directory', default = None)

	args = parser.parse_args()

	_size = None

	_size = read_size(args.sizeDir)
	edit_yaml(args.yamlDir, _size)







