#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(
    description="single cell RNAseq mixed genotype clustering using sparse mixture model clustering with tensorflow.")
parser.add_argument("-i","--bam",required=True, help="cellranger bam")
parser.add_argument("-b","--barcodes",required=True, help= "barcodes.tsv from cellranger")
parser.add_argument("-f","--fasta",required=True, help="reference fasta file")
parser.add_argument("-t","--threads",required=True, help="max threads to use")
parser.add_argument("--ignore",required=False, default="False", help = "set to True to ignore data error assertions")
args = parser.parse_args()

# importing all reqs to make sure things are installed
import numpy as np
import tensorflow as tf
import scipy
import math
import pystan
import pyvcf
import pyfasta
import subprocess

#load each file to make sure it is legit
bc_set = set()
with open(args.barcodes) as barcodes:
    for (index, line) in enumerate(barcodes):
        bc = line.strip()
        bc_set.add(bc)

assert(len(bc_set) > 50, "Fewer than 50 barcodes in barcodes file? We expect 1 barcode per line.")

#test bam load
bam = pysam.AlignmentFile(args.bam)
num_cb = 0
num_cb_cb = 0 # num reads with barcodes from barcodes.tsv file
num_umi = 0
num_read_test = 100000
for (index,read) in enumerate(bam):
    if index >= num_read_test:
        break
    if read.has_tag("CB"):
        num_cb += 1
        if read.get_tag("CB") in bc_set:
            num_cb_cb += 1
    if read.has_tag("UB"):
        num_umi += 1
if not args.ignore == "True":
    assert(float(num_cb)/float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have cell barcode tag (CB), turn on --ignore True to ignore")
    assert(float(num_umi)/float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have UMI tag (UB), turn on --ignore True to ignore")
    assert(float(num_cb_cb)/float(num_read_test) > 0.25, "Less than 25% of first 100000 reads have cell barcodes from barcodes file, is this the correct barcode file? turn on --ignore True to ignore"_
#test fasta load
fasta = pyfasta.Fasta(args.fasta)

bam = pysam.AlignmentFile(args.bam)
total_reference_length = 0
for chrom in bam.references:
    total_reference_length += bam.get_reference_length(chrom)
step_length = total_reference_length/int(args.threads)
regions = []
region = []
region_so_far = 0
chrom_so_far = 0
for chrom in bam.references:
    chrom_length = bam.get_reference_length(chrom)
    while True:
        if region_so_far + chrom_length < step_length:
            region.append((chrom, chrom_so_far, chrom_length))
            region_so_far += chrom_length - chrom_so_far
            chrom_so_far = 0
            break
        else:
            region.append((chrom, chrom_so_far, step_length - region_so_far))
            regions.append(region)
            region = []
            chrom_so_far = step_length - region_so_far + 1
            region_so_far = 0
if len(region) > 0:
    regions.append(region)

# for testing, delete this later
print(len(regions))
print(regions)

procs = []
for region in regions:
    
