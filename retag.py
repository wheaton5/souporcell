#!/usr/bin/env python

import pysam
import argparse


parser = argparse.ArgumentParser(description = 'Retag reads with their cell barcodes and UMIs')

parser.add_argument('-s', '--sam', required = True, help = "sam file")
parser.add_argument('-o', '--out', required = True, help = "output (will be a bam)")
parser.add_argument("--no_umi", required = False, default = False, type = bool, help = "set True if your bam has no umi tag")
parser.add_argument("--umi_tag", required = False, default = "UB", help = "set if umi tag in output bam should not be UB")
parser.add_argument("--cell_tag", required = False, default = "CB", help = "set if cell barcode tag should not be CB")
args = parser.parse_args()

CELL_TAG = args.cell_tag
UMI_TAG = args.umi_tag

bam = pysam.AlignmentFile(args.sam)

bamout = pysam.AlignmentFile(args.out,'wb', template = bam)

for read in bam:
    qname = read.qname
    tokens = qname.split(";")
    if args.no_umi:
        assert len(tokens) == 2
        read.set_tag(CELL_TAG, tokens[-1])
    else:
        assert len(tokens) == 3
        read.set_tag(CELL_TAG, tokens[-2])
        read.set_tag(UMI_TAG, tokens[-1])
    bamout.write(read)

bamout.close()
