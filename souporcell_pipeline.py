#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(
    description="single cell RNAseq mixed genotype clustering using sparse mixture model clustering with tensorflow.")
parser.add_argument("-i","--bam",required=True, help="cellranger bam")
parser.add_argument("-b","--barcodes",required=True, help= "barcodes.tsv from cellranger")
parser.add_argument("-f","--fasta",required=True, help="reference fasta file")
parser.add_argument("-t","--threads",required=True, help="max threads to use")
parser.add_argument("-o","--out_dir",required=True, help="name of directory to place souporcell files")
parser.add_argument("-k","--clusters",required=True, help="number cluster, tbd add easy way to run on a range of k")
parser.add_argument("-p","--ploidy",required=False,default="2",help="ploidy, must be 1 or 2, default = 2")
parser.add_argument("--min_alt",required=False, default="10", help="min alt to use locus, default = 10.")
parser.add_argument("--min_ref",required=False, default="10", help="min ref to use locus, default = 10.")
parser.add_argument("--max_loci",required=False, default="2048",help="max loci per cell, affects speed, default = 2048.")
parser.add_argument("--ignore",required=False, default="False", help = "set to True to ignore data error assertions")
args = parser.parse_args()


print("checking modules")
# importing all reqs to make sure things are installed
import numpy as np
import tensorflow as tf
import scipy
import gzip
import math
import pystan
import vcf
import pysam
import pyfasta
import subprocess
import time
import os
print("imports done")
subprocess.check_call(["mkdir",args.out_dir])

print("checking bam for expected tags")
#load each file to make sure it is legit
bc_set = set()
with open(args.barcodes) as barcodes:
    for (index, line) in enumerate(barcodes):
        bc = line.strip()
        bc_set.add(bc)

assert len(bc_set) > 50, "Fewer than 50 barcodes in barcodes file? We expect 1 barcode per line."

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
    assert float(num_cb)/float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have cell barcode tag (CB), turn on --ignore True to ignore"
    assert float(num_umi)/float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have UMI tag (UB), turn on --ignore True to ignore"
    assert float(num_cb_cb)/float(num_read_test) > 0.05, "Less than 25% of first 100000 reads have cell barcodes from barcodes file, is this the correct barcode file? turn on --ignore True to ignore"
print("checking fasta")
#test fasta load
fasta = pyfasta.Fasta(args.fasta)

bam = pysam.AlignmentFile(args.bam)
total_reference_length = 0
for chrom in bam.references:
    total_reference_length += bam.get_reference_length(chrom)
step_length = int(math.ceil(total_reference_length/int(args.threads)))
regions = []
region = []
region_so_far = 0
chrom_so_far = 0
print("creating chunks")
for chrom in bam.references:
    chrom_length = bam.get_reference_length(chrom)
    while True:
        if region_so_far + (chrom_length - chrom_so_far) < step_length:
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
args.threads = int(args.threads)
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
                assert not(procs[index].returncode),"renamer subprocess terminated abnormally with code "+str(procs[index].returncode)
        if len(region_fastqs[index]) == len(region) - 1:
            block = True
        if not block:
            sub_index = len(region_fastqs[index])
            chrom = region[sub_index][0]
            start = region[sub_index][1]
            end = region[sub_index][2]
            fq_name = args.out_dir+"/souporcell_fastq_"+str(index)+"_"+str(sub_index)+".fq"
            p = subprocess.Popen(["renamer.py","--bam",args.bam,"--barcodes",args.barcodes,"--out",fq_name,"--chrom",chrom,"--start",str(start),"--end",str(end)])
            all_fastqs.append(fq_name)
            procs[index] = p
            region_fastqs[index].append(fq_name)
            any_running = True
    time.sleep(20)
    print("tick "+str(any_running)) 
 
print("remapping with minimap2")
FNULL = open(os.devnull, 'w') # minimap too verbose
# run minimap2
minimap_tmp_files = []
for index in range(args.threads):
    output = args.out_dir+"/souporcell_minimap_tmp_"+str(index)+".sam"
    minimap_tmp_files.append(output)
    with open(output, 'w') as samfile:
        with open(args.out_dir+"/minimap.err",'w') as minierr:
            subprocess.check_call(["minimap2","-ax","splice","-t",str(args.threads),"-G50k","-k","21","-w","11","--sr","-A2","-B8","-O12,32","-E2,1",
                "-r200","-p.5","-N20","-f1000,5000","-n2","-m20","-s40","-g2000","-2K50m","--secondary=no",args.fasta]+region_fastqs[index], stdout = samfile, stderr = minierr)

print("cleaning up tmp fastqs")
# clean up tmp fastqs
for fq in all_fastqs:
    subprocess.check_call(["rm",fq])

print("repopulating cell barcode and UMI tags")
# run retagger
procs = []
retag_files = []
for index in range(args.threads):
    outfile = args.out_dir+"/souporcell_retag_tmp_"+str(index)+".bam"
    retag_files.append(outfile)
    p = subprocess.Popen(["retag.py","--sam",minimap_tmp_files[index],"--out",outfile])
    procs.append(p)
for p in procs: # wait for processes to finish
    p.wait()
    assert not(p.returncode),"retag subprocess ended abnormally with code "+str(p.returncode)

print("cleaning up tmp samfiles")
# clean up tmp samfiles
for samfile in minimap_tmp_files:
    subprocess.check_call(["rm",samfile])

print("sorting retagged bam files")
# sort retagged files
sort_jobs = []
file_handles = []
filenames = []
for index in range(args.threads):
    filename = args.out_dir+"/souporcell_retag_sorted_tmp_"+str(index)+".bam"
    filenames.append(filename)
    filehandle = open(filename,'wb')
    file_handles.append(filehandle)
    p = subprocess.Popen(["samtools","sort",retag_files[index]],stdout=filehandle,stderr=FNULL)
    sort_jobs.append(p)
    
# wait for jobs to finish
for job in sort_jobs:
    job.wait()
    assert not(job.returncode),"samtools sort ended abnormally with code "+str(job.returncode)
#close files
for filehandle in file_handles:
    filehandle.close()

#clean up unsorted bams
for bam in retag_files:
    subprocess.check_call(["rm",bam])

print("merging sorted bams")
final_bam = args.out_dir+"/souporcell_minimap_tagged_sorted.bam"
subprocess.check_call(["samtools","merge", final_bam]+filenames)

subprocess.check_call(["samtools","index", final_bam])

# clean up tmp bams
for filename in filenames:
    subprocess.check_call(['rm',filename])


total_reference_length = 0
for chrom in sorted(fasta.keys()):
    total_reference_length += len(fasta[chrom])

regions = []
region = []
region_so_far = 0
chrom_so_far = 0
for chrom in sorted(fasta.keys()):
    chrom_length = len(fasta[chrom])
    while True:
        if region_so_far + (chrom_length - chrom_so_far) < step_length:
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

region_vcfs = [[] for x in range(args.threads)]
all_vcfs = []
procs = [None for x in range(args.threads)]
any_running = True
filehandles = []
# run renamer in parallel manner
print("running freebayes")
while any_running:
    any_running = False
    for (index, region) in enumerate(regions):
        block = False
        if procs[index]:
            block = procs[index].poll() == None
            if block:
                any_running = True
            else:
                assert not(procs[index].returncode),"freebayes subprocess terminated abnormally with code "+str(procs[index].returncode)
        if len(region_vcfs[index]) == len(region) - 1:
            block = True
        if not block:
            sub_index = len(region_vcfs[index])
            chrom = region[sub_index][0]
            start = region[sub_index][1]
            end = region[sub_index][2]
            vcf_name = args.out_dir+"/souporcell_"+str(index)+"_"+str(sub_index)+".vcf"
            filehandle = open(vcf_name,'w')
            filehandles.append(filehandle)
            p = subprocess.Popen(["freebayes","-f",args.fasta,"-r",chrom+":"+str(start)+"-"+str(end),"-iXu","-C","2",
                "-q","20","-n", "3", "-E","1","-m","30","--min-coverage","6","--max-coverage","100000","--pooled-continuous",final_bam],stdout=filehandle)
            all_vcfs.append(vcf_name)
            procs[index] = p
            region_vcfs[index].append(vcf_name)
            any_running = True
    time.sleep(10)
for filehandle in filehandles:
    filehandle.close()
print("merging vcfs")
with open(args.out_dir+"/souporcell_merged_vcf.vcf",'w') as vcfout:
    subprocess.check_call(["bcftools","concat"]+all_vcfs,stdout=vcfout)
subprocess.check_call(['rm']+all_vcfs)
with open(args.out_dir+"/souporcell_merged_sorted_vcf.vcf",'w') as vcfout:
    subprocess.check_call(['bcftools','sort',args.out_dir+"/souporcell_merged_vcf.vcf"],stdout=vcfout,stderr=FNULL)
subprocess.check_call(['rm',args.out_dir+'/souporcell_merged_vcf.vcf'])
subprocess.check_call(['bgzip',args.out_dir+"/souporcell_merged_sorted_vcf.vcf"])
final_vcf = args.out_dir+"/souporcell_merged_sorted_vcf.vcf.gz"
subprocess.check_call(['tabix','-p','vcf',final_vcf])

print("running vartrix")
ref_mtx = args.out_dir+"/ref.mtx"
alt_mtx = args.out_dir+"/alt.mtx"
subprocess.check_call(["vartrix","--umi","--mapq","30","-b",final_bam,"-c",args.barcodes,"--scoring-method","coverage","--threads",str(args.threads),
    "--ref-matrix",ref_mtx, "--out-matrix",alt_mtx,"-v",final_vcf, "--fasta", args.fasta])

print("running souporcell clustering")
cluster_file = args.out_dir+"/clusters_tmp.tsv"
with open(args.out_dir+"/souporcell.log",'w') as log:
    subprocess.check_call(["souporcell.py","-a",alt_mtx,"-r",ref_mtx,"-b",args.barcodes,"-k",args.clusters,
        "-t",str(args.threads), "-l", args.max_loci, "--min_alt", args.min_alt, "--min_ref", args.min_ref,'--out',cluster_file],stdout=log,stderr=log) 

print("running souporcell doublet detection")
doublet_file = args.out_dir + "/clusters.tsv"
with open(doublet_file,'w') as dub:
    subprocess.check_call(["troublet","--alts",alt_mtx,"--refs",ref_mtx,"--clusters",cluster_file],stdout=dub)

print("running co inference of ambient RNA and cluster genotypes")
subprocess.check_call(["consensus.py","-c",doublet_file,"-a",alt_mtx,"-r",ref_mtx,"-p", args.ploidy,
    "--soup_out",args.out_dir+"/ambient_rna.txt","--vcf_out",args.out_dir+"/cluster_genotypes.vcf","--vcf",final_vcf])
