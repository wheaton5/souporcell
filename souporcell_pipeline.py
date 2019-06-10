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
import time

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
region_fastqs = [[] for x in range(args.threads)]
all_fastqs = []
procs = [None for x in range(args.threads)]
any_running = True
# run renamer in parallel manner
print("generating fastqs with cell barcodes and umis in readname")
while any_running:
    any_running = False
    for (index, region) in enumerate(regions):
        block = False
        if procs[index]:
            block = procs[index].poll() == None
            if block:
                any_running = True
            else:
                assert(not(procs[index].returncode),"renamer subprocess terminated abnormally with code "+str(procs[index].returncode))
        if len(region_fastqs) == len(region):
            block = True
        if not block:
            sub_index = len(region_fastqs[index])
            chrom = region[index][sub_index][0]
            start = region[index][sub_index][1]
            end = region[index][sub_index][2]
            fq_name = "souporcell_fastq_"+str(index)+"_"+str(sub_index)+".fq"
            p = subprocess.Popen(["renamer.py","--bam",args.bam,"--barcodes",args.barcodes,"--out",fq_name,"--chrom",chrom,"--start",str(start),"--end",str(end)])
            all_fastqs.append(fq_name)
            procs[index] = p
            region_fastqs[index].append(fq_name)
            any_running = True
    time.sleep(20)

print("remapping with minimap2")
# run minimap2
minimap_tmp_files = []
for index in range(args.threads):
    output = "souporcell_minimap_tmp_"+str(index)+".sam")
    minimap_tmp_files.append(output)
    with open(output, 'w') as samfile:
        subprocess.check_call(["minimap2","-ax","splice","-t",str(args.threads),"-G50k","-k","21","-w","11","--sr","-A2","-B8","-O12,32","-E2,1",
            "-r200","-p.5","-N20","-f1000,5000","-n2","-m20","-s40","-g2000","-2K50m","--secondary=no",args.fasta]+region_fastqs, stdout = samfile)

print("cleaning up tmp fastqs")
# clean up tmp fastqs
for fq in all_fastqs:
    suprocess.check_call(["rm",fq])

print("repopulating cell barcode and UMI tags")
# run retagger
procs = []
retag_files = []
for index in range(args.threads):
    outfile = "souporcell_retag_tmp_"+str(index)+".bam"
    retag_files.append(outfile)
    p = subprocess.Popen(["retag.py","--sam",minimap_tmp_files[index],"--out",outfile)
    procs.append(p)
for p in procs: # wait for processes to finish
    p.join()
    assert(not(p.returncode),"retag subprocess ended abnormally with code "+str(p.returncode))

print("cleaning up tmp samfiles")
# clean up tmp samfiles
for samfile in minimap_tmp_files:
    subprocess.check_call(["rm",samfile])

print("sorting retagged bam files")
# sort retagged files
for index in range(args.threads):
    


