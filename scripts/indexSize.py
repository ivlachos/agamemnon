#!/usr/bin/env python3

import argparse
import os



def calculate_size(indexDir):
	size = 0
	for dirpath, dirnames, filenames in os.walk(indexDir):
		for f in filenames:
			fp = os.path.join(dirpath, f)
			size += os.path.getsize(fp)
	size = int(size / (10 ** 6)) + 3

	return size



def write_size(size, outDir):
	results = open(outDir + '/indexSize.tmp', 'w')
	results.write(str(size))

	return None



if __name__ == '__main__':
	parser = argparse.ArgumentParser(__file__, description = "Calculate Index Size")
	parser.add_argument('-indexDir', '-id', help = 'Index directory', default = None)
	parser.add_argument('-outDir', '-od', help = 'Output directory', default = None)

	args = parser.parse_args()

	_size = 0

	_size = calculate_size(args.indexDir)
	write_size(_size, args.outDir)



