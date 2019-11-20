#!/usr/bin/env python

import pysam
import argparse

parser = argparse.ArgumentParser(description='make fastq from possorted_genome_bam.bam from cellranger')

parser.add_argument('-f', '--bam', required=True, help="cellranger bam")
parser.add_argument('-b', '--barcodes', required=True, help="cellranger barcodes.tsv")
parser.add_argument('-o', '--out', required=True, help="output fastq name")
parser.add_argument('-c', '--chrom', required = False, help="chrom")
parser.add_argument('-s', '--start', required = False, help="start")
parser.add_argument('-e', '--end', required = False, help="end")
args = parser.parse_args()

assert (not(args.chrom) and not(args.start) and not(args.end)) or (args.chrom and args.start and args.end), "if specifying region, must specify chrom, start, and end"

fn = args.bam#"possorted_genome_bam.bam"#files[0]
bam = pysam.AlignmentFile(fn, "rb")

cell_barcodes = set([])
with open(args.barcodes) as barcodes:
    for line in barcodes:
        tokens=line.strip().split()
        cell_barcodes.add(tokens[0])

if args.chrom:
    bam = bam.fetch(args.chrom, int(args.start), int(args.end))

recent_umis = {}
with open(args.out,'w') as fastq:
    for (index,read) in enumerate(bam):
        if not read.has_tag("CB"):
            continue
        cell_barcode = read.get_tag("CB")
        if read.is_secondary or read.is_supplementary:
            continue
        if not read.has_tag("UB"):
            continue
        UMI = read.get_tag("UB")
        pos = read.pos
        full_umi = cell_barcode + UMI + str(pos)
        if full_umi in recent_umis:
            continue
        #recent_umis[full_umi] = pos
        #if index % 10000 == 0:
        #    keys_to_remove = []
        #    for (key, val) in recent_umis.items():
        #        if val - pos > 200:
        #            keys_to_remove.append(key)
        #    for key in keys_to_remove:
        #        del recent_umis[key]
    
        readname = read.qname
        if read.has_tag("CB") and read.get_tag("CB") in cell_barcodes:
            fastq.write("@"+read.qname+";"+cell_barcode+";"+UMI+"\n")
            fastq.write(read.seq+"\n")
            fastq.write("+\n")
            fastq.write(read.qual+"\n")

