#!/usr/bin/env python

##Written by Haynes Heaton to reformat SS2 files to mimic 10x format for Souporcell


import argparse

parser = argparse.ArgumentParser(description = "convert smartseq data to a format souporcell can use")
parser.add_argument("-i", "--input_dir", required = True, help = "directory with smartseq bams")
parser.add_argument("-o", "--output_prefix", required= True, help = "outputs bam and barcode files with this prefix")
args = parser.parse_args()

import glob

import pysam

bamfiles = glob.glob(args.input_dir+"/*.BAM")
assert len(bamfiles) > 0, "I don't see bam files in directory" + args.input_dir
template = pysam.AlignmentFile(bamfiles[0])
outbam = pysam.AlignmentFile(args.output_prefix+".bam",'wb', template = template)
cell_barcode = 0

with open(args.output_prefix +".tsv", "w") as barcodes:
    for (cell_barcode, bamfile) in enumerate(bamfiles):
        umi = 0
        inbam = pysam.AlignmentFile(bamfile)
        barcode = bamfile.split('/')[-1]
        barcodes.write(barcode+"\n")
        for read in inbam:
            read.set_tag("UB", str(umi))
            umi += 1
            read.set_tag("CB", str(barcode))
            outbam.write(read)

import subprocess

with open(args.output_prefix+"_sorted.bam",'w') as bamout2:
    subprocess.check_call(["samtools","sort",args.output_prefix+".bam"],stdout = bamout2)

subprocess.check_call(["rm", args.output_prefix+".bam"])
subprocess.check_call(["mv", args.output_prefix+"_sorted.bam",args.output_prefix+".bam"])

