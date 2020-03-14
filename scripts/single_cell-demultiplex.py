#!/usr/bin/env python3

from gatb import Bank
import argparse


def barcodes_dict(bindex):
    tmp_barcodes = {}
    for seq in bindex:
        tmp_barcodes[seq.comment.decode('utf-8').split(' ')[0].split('-')[1].strip()] = \
            [seq.sequence.decode('utf-8').strip(), seq.quality.decode('utf-8')]

    return tmp_barcodes



def groupReads(bm, barcodes, fileName):
    barcodesReads = {}

    for seq in bm:
        if seq.comment.decode('utf-8').split(' ')[0].split('-')[1].strip() in barcodes:
            val = qualityControl(barcodes[seq.comment.decode('utf-8').split(' ')[0].split('-')[1].strip()][1])
            if val:
                if barcodes[seq.comment.decode('utf-8').split(' ')[0].split('-')[1].strip()][0] in barcodesReads:
                    barcodesReads[barcodes[seq.comment.decode('utf-8').split(' ')[0].split('-')[1].strip()][0]].\
                        append([seq.comment.decode('utf-8'), seq.sequence.decode('utf-8'), '+', seq.quality.decode('utf-8')])
                else:
                    barcodesReads[barcodes[seq.comment.decode('utf-8').split(' ')[0].split('-')[1].strip()][0]] = \
                        [[seq.comment.decode('utf-8'), seq.sequence.decode('utf-8'), '+', seq.quality.decode('utf-8')]]
                fileName.write('@' + seq.comment.decode('utf-8') + '-' \
			+ barcodes[seq.comment.decode('utf-8').split(' ')[0].split('-')[1].strip()][0] + '\n')
                fileName.write(seq.sequence.decode('utf-8') + '\n')
                fileName.write('+\n')
                fileName.write(seq.quality.decode('utf-8') + '\n')

    return barcodesReads



def qualityControl(read_quality):
    passVal = True
    for char in read_quality:
        if Qscore[char] < barQ:
            passVal = False
            break

    return passVal



def writeFastqFiles(barcodesReads):
    num = 0
    for key, value in barcodesReads.items():
        if len(value) >= int(args.threshold):
            current_file = open(args.out + key + '.fq', 'a+')
            for i in range(0, len(value)):
                num += 1
                current_file.write('@' + key + '-' + str(num) + '-' + value[i][0] + '\n')
                current_file.write(value[i][1] + '\n' + value[i][2] + '\n' + value[i][3] + '\n')
            current_file.close()
        else:
            continue

    return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__file__, description='Single-cell fastq files handler')
    parser.add_argument('-m1Fastq', '-m1', help = 'Fastq file mate 1', default = None)
    parser.add_argument('-m2Fastq', '-m2', help = 'Fastq file mate 2 (if paired-end)', default = None)
    parser.add_argument('-paired', '-pe', help = 'paired-end [0 - 1]', default = 0)
    parser.add_argument('-iFastq', '-if', help = 'Indexed Fastq file', default = None)
    parser.add_argument('-mode', '-mo', help = '''Mode 1 [Create one Fastq file passing barcodes in readnames] | 
			Mode 2 [Create multiple Fastq files, one for each barcode]''', default = 1)
    parser.add_argument('-threshold', '-th', help = 'Single-cell cutoff size', default = 50)
    parser.add_argument('-out', '-o', help = 'Output directory', default = None)
    args = parser.parse_args()

    bank_index = Bank(args.iFastq)
    bank_m1 = Bank(args.m1Fastq)
    if int(args.paired) == 1:
        bank_m2 = Bank(args.m2Fastq)
    barcodes = {}
    barcodesReads = {}
    barQ = 15
    Qscore = dict((chr(i), i - 33) for i in range(33, 90))
    if int(args.paired) == 0:
        fqFile = open(args.out + 'barcodedReads.fq', 'w')
    else:
        fqFile1 = open(args.out + 'barcodedReads_1.fq', 'w')
        fqFile2 = open(args.out + 'barcodedReads_2.fq', 'w')

    barcodes = barcodes_dict(bank_index)
    if int(args.paired) == 0:
        barcodesReads = groupReads(bank_m1, barcodes, fqFile)
    elif int(args.paired) == 1:
        barcodesReads = groupReads(bank_m1, barcodes, fqFile1)
        barcodesReads = groupReads(bank_m2, barcodes, fqFile2)

    if int(args.mode) == 2:
        writeFastqFiles(barcodesReads)


