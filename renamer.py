#!/usr/bin/env python

import pysam
import argparse
import gzip
parser = argparse.ArgumentParser(description='make fastq from possorted_genome_bam.bam from cellranger')

parser.add_argument('-f', '--bam', required=True, help="cellranger bam")
parser.add_argument('-b', '--barcodes', required=True, help="cellranger barcodes.tsv")
parser.add_argument('-o', '--out', required=True, help="output fastq name")
parser.add_argument('-c', '--chrom', required = False, help="chrom")
parser.add_argument('-s', '--start', required = False, help="start")
parser.add_argument('-e', '--end', required = False, help="end")
parser.add_argument("--no_umi", required = False, default = "False", help = "set to True if your bam has no umi tag")
parser.add_argument("--umi_tag", required = False, default = "UB", help = "set if umi tag is not UB")
parser.add_argument("--cell_tag", required = False, default = "CB", help = "set if cell barcode tag is not CB")
args = parser.parse_args()

if args.no_umi == "True":
    args.no_umi = True
else:
    args.no_umi = False


UMI_TAG = args.umi_tag
CELL_TAG = args.cell_tag

assert (not(args.chrom) and not(args.start) and not(args.end)) or (args.chrom and args.start and args.end), "if specifying region, must specify chrom, start, and end"

fn = args.bam#"possorted_genome_bam.bam"#files[0]
bam = pysam.AlignmentFile(fn, "rb")

open_function = lambda f: gzip.open(f,"rt") if f[-3:] == ".gz" else open(f)
cell_barcodes = set([])
with open_function(args.barcodes) as barcodes:
    for line in barcodes:
        tokens=line.strip().split()
        cell_barcodes.add(tokens[0])

if args.chrom:
    bam = bam.fetch(args.chrom, int(args.start), int(args.end))

recent_umis = {}
with open(args.out,'w') as fastq:
    for (index,read) in enumerate(bam):
        if not read.has_tag(CELL_TAG):
            continue
        cell_barcode = read.get_tag(CELL_TAG)
        if read.is_secondary or read.is_supplementary:
            continue
        pos = read.pos
        if args.no_umi:
            full_umi = cell_barcode + str(pos)
        else:
            if not read.has_tag(UMI_TAG):
                continue
            UMI = read.get_tag(UMI_TAG)
            full_umi = cell_barcode + UMI + str(pos)
        if full_umi in recent_umis:
            continue
        if read.seq is None:
            continue
    
        readname = read.qname
        if read.has_tag(CELL_TAG) and read.get_tag(CELL_TAG) in cell_barcodes:
            if args.no_umi:
                fastq.write("@"+read.qname+";"+cell_barcode+"\n")
            else:
                fastq.write("@"+read.qname+";"+cell_barcode+";"+UMI+"\n")
            fastq.write(read.seq+"\n")
            fastq.write("+\n")
            fastq.write(read.qual+"\n")

