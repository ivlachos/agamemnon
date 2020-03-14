#! /usr/bin/env python3
"""
Take a file full of accessions and one or more NCBI 'accession2taxid' files
and create a CSV named 'acc_file.taxid' that contains the accessions
and their associated taxids.
"""

from __future__ import print_function
import argparse
import gzip

def main():
    p = argparse.ArgumentParser()
    p.add_argument('outDir')
    p.add_argument('acc_file')
    p.add_argument('acc2taxid_files', nargs='+')
    args = p.parse_args()

    with open(args.acc_file) as fp:
        acc_set = set([ x.strip() for x in fp ])

    outfp = open(args.outDir + 'accessionsTaxIDs.tab', 'w')

    m = 0
    for filename in args.acc2taxid_files:
        if not acc_set: break

        xopen = open
        if filename.endswith('.gz'):
            xopen = gzip.open

        with xopen(filename, 'rt') as fp:
            next(fp)                #  skip headers
            for n, line in enumerate(fp):
                if not acc_set: break

                try:
                    _, acc_version, taxid, _ = line.split()
                except ValueError:
                    continue

                if acc_version in acc_set:
                    m += 1
                    outfp.write('{}\t{}\n'.format(acc_version, taxid))
                    acc_set.remove(acc_version)


if __name__ == '__main__':
    main()
