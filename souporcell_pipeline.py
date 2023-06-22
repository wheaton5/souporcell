#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(
    description="single cell RNAseq mixed genotype clustering using sparse mixture model clustering.")
parser.add_argument("-i", "--bam", required = True, help = "cellranger bam")
parser.add_argument("-b", "--barcodes", required = True, help = "barcodes.tsv from cellranger")
parser.add_argument("-f", "--fasta", required = True, help = "reference fasta file")
parser.add_argument("-t", "--threads", required = True, type = int, help = "max threads to use")
parser.add_argument("-o", "--out_dir", required = True, help = "name of directory to place souporcell files")
parser.add_argument("-k", "--clusters", required = True, help = "number cluster, tbd add easy way to run on a range of k")
parser.add_argument("-p", "--ploidy", required = False, default = "2", help = "ploidy, must be 1 or 2, default = 2")
parser.add_argument("--min_alt", required = False, default = "10", help = "min alt to use locus, default = 10.")
parser.add_argument("--min_ref", required = False, default = "10", help = "min ref to use locus, default = 10.")
parser.add_argument("--max_loci", required = False, default = "2048", help = "max loci per cell, affects speed, default = 2048.")
parser.add_argument("--restarts", required = False, default = 100, type = int, 
    help = "number of restarts in clustering, when there are > 12 clusters we recommend increasing this to avoid local minima")
parser.add_argument("--common_variants", required = False, default = None, 
    help = "common variant loci or known variant loci vcf, must be vs same reference fasta")
parser.add_argument("--known_genotypes", required = False, default = None, 
    help = "known variants per clone in population vcf mode, must be .vcf right now we dont accept gzip or bcf sorry")
parser.add_argument("--known_genotypes_sample_names", required = False, nargs = '+', default = None, 
    help = "which samples in population vcf from known genotypes option represent the donors in your sample")
parser.add_argument("--skip_remap", required = False, default = False, type = bool, 
    help = "don't remap with minimap2 (not recommended unless in conjunction with --common_variants")
parser.add_argument("--no_umi", required = False, default = "False", help = "set to True if your bam has no UMI tag, will ignore/override --umi_tag")
parser.add_argument("--umi_tag", required = False, default = "UB", help = "set if your umi tag is not UB")
parser.add_argument("--cell_tag", required = False, default = "CB", help = "DOES NOT WORK, vartrix doesnt support this! set if your cell barcode tag is not CB")
parser.add_argument("--ignore", required = False, default = False, type = bool, help = "set to True to ignore data error assertions")
parser.add_argument("--aligner", required = False, default = "minimap2", help = "optionally change to HISAT2 if you have it installed, not included in singularity build")
args = parser.parse_args()



if args.no_umi == "True":
    args.no_umi = True
else:
    args.no_umi = False

print("checking modules")
# importing all reqs to make sure things are installed
import numpy as np
import scipy
import gzip
import math
import pystan
import vcf
import pysam
import pyfaidx
import subprocess
import time
import os
print("imports done")


open_function = lambda f: gzip.open(f,"rt") if f[-3:] == ".gz" else open(f)

print("checking bam for expected tags")
UMI_TAG = args.umi_tag
CELL_TAG = args.cell_tag
assert CELL_TAG == "CB", "vartrix doesnt support different cell tags, remake bam with cell tag as CB"
#load each file to make sure it is legit
bc_set = set()
with open_function(args.barcodes) as barcodes:
    for (index, line) in enumerate(barcodes):
        bc = line.strip()
        bc_set.add(bc)

assert len(bc_set) > 50, "Fewer than 50 barcodes in barcodes file? We expect 1 barcode per line."
assert args.aligner == "minimap2" or args.aligner == "HISAT2", "--aligner expects minimap2 or HISAT2"

assert not(not(args.known_genotypes == None) and not(args.common_variants == None)), "cannot set both know_genotypes and common_variants"
if args.known_genotypes_sample_names:
    assert not(args.known_genotypes == None), "if you specify known_genotype_sample_names, must specify known_genotypes option"
    assert len(args.known_genotypes_sample_names) == int(args.clusters), "length of known genotype sample names should be equal to k/clusters"
if args.known_genotypes:
    reader = vcf.Reader(open(args.known_genotypes))
    assert len(reader.samples) >= int(args.clusters), "number of samples in known genotype vcfs is less than k/clusters"
    if args.known_genotypes_sample_names == None:
        args.known_genotypes_sample_names = reader.samples
    for sample in args.known_genotypes_sample_names:
        assert sample in args.known_genotypes_sample_names, "not all samples in known genotype sample names option are in the known genotype samples vcf?"

#test bam load
bam = pysam.AlignmentFile(args.bam)
num_cb = 0
num_cb_cb = 0 # num reads with barcodes from barcodes.tsv file
num_umi = 0
num_read_test = 100000
for (index,read) in enumerate(bam):
    if index >= num_read_test:
        break
    if read.has_tag(CELL_TAG):
        num_cb += 1
        if read.get_tag(CELL_TAG) in bc_set:
            num_cb_cb += 1
    if read.has_tag(UMI_TAG):
        num_umi += 1
if not args.ignore:
    if args.skip_remap and args.common_variants == None and args.known_genotypes == None:
        assert False, "WARNING: skip_remap enables without common_variants or known genotypes. Variant calls will be of poorer quality. Turn on --ignore True to ignore this warning"
        
    assert float(num_cb) / float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have cell barcode tag (CB), turn on --ignore True to ignore"
    if not(args.no_umi):
        assert float(num_umi) / float(num_read_test) > 0.5, "Less than 50% of first 100000 reads have UMI tag (UB), turn on --ignore True to ignore"
    assert float(num_cb_cb) / float(num_read_test) > 0.05, "Less than 25% of first 100000 reads have cell barcodes from barcodes file, is this the correct barcode file? turn on --ignore True to ignore"

print("checking fasta")
fasta = pyfaidx.Fasta(args.fasta, key_function = lambda key: key.split()[0])

def get_fasta_regions(fastaname, threads):
    fasta = pyfaidx.Fasta(args.fasta, key_function = lambda key: key.split()[0])
    total_reference_length = 0
    for chrom in sorted(fasta.keys()):
       total_reference_length += len(fasta[chrom])
    step_length = int(math.ceil(total_reference_length/threads))
    regions = []
    region = []
    region_so_far = 0
    chrom_so_far = 0
    for chrom in sorted(fasta.keys()):
        chrom_length = len(fasta[chrom])
        chrom_so_far = 0
        if chrom_length < 250000:
            continue
        while True:
            if region_so_far + (chrom_length - chrom_so_far) < step_length:
                region.append((chrom, chrom_so_far, chrom_length))
                region_so_far += chrom_length - chrom_so_far
                chrom_so_far = 0
                break
            else:
                region.append((chrom, chrom_so_far, chrom_so_far + step_length - region_so_far))
                regions.append(region)
                region = []
                chrom_so_far += step_length - region_so_far
                region_so_far = 0
    if len(region) > 0:
        if len(regions) == args.threads:
            regions[-1] = regions[-1] + region
        else:
            regions.append(region)
    return regions


def get_bam_regions(bamname, threads):
    bam = pysam.AlignmentFile(bamname)
    total_reference_length = 0
    for chrom in bam.references:
        total_reference_length += bam.get_reference_length(chrom)
    step_length = int(math.ceil(total_reference_length / threads))
    regions = []
    region = []
    region_so_far = 0
    chrom_so_far = 0
    for chrom in bam.references:
        chrom_length = bam.get_reference_length(chrom)
        #print(chrom+" size "+str(chrom_length)+" and step size "+str(step_length))
        while True:
            #print("\tregion so far "+str(region_so_far)+" chrom so far "+str(chrom_so_far)) 
            #print("\t"+str(chrom_length - chrom_so_far)+" <= "+str(step_length - region_so_far))
            #print("\t"+str((chrom_length - chrom_so_far) <= step_length - region_so_far))
            if (chrom_length - chrom_so_far) <= step_length - region_so_far:
                region.append((chrom, chrom_so_far, chrom_length))
                #print("\t\tending chrom\t"+chrom+":"+str(chrom_so_far)+"-"+str(chrom_length))
                region_so_far += chrom_length - chrom_so_far
                chrom_so_far = 0
                break
            else:
                region.append((chrom, chrom_so_far, chrom_so_far + step_length - region_so_far))
                #print("\t\tending region\t"+chrom+":"+str(chrom_so_far)+"-"+str(chrom_so_far + step_length - region_so_far))
                regions.append(region)
                region = []
                chrom_so_far += step_length - region_so_far
                region_so_far = 0
    if len(region) > 0:
        regions.append(region)
    
    return regions

def make_fastqs(args):
    if not os.path.isfile(args.bam + ".bai"):
        print("no bam index found, creating")
        subprocess.check_call(['samtools', 'index', args.bam])
    if not os.path.isfile(args.fasta + ".fai"):
        print("fasta index not found, creating")
        subprocess.check_call(['samtools', 'faidx', args.fasta])
    regions = get_bam_regions(args.bam, int(args.threads))
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
                    assert not(procs[index].returncode), "renamer subprocess terminated abnormally with code " + str(procs[index].returncode)
            if len(region_fastqs[index]) == len(region):
                block = True
            if not block:
                sub_index = len(region_fastqs[index])
                chrom = region[sub_index][0]
                start = region[sub_index][1]
                end = region[sub_index][2]
                fq_name = args.out_dir + "/souporcell_fastq_" + str(index) + "_" + str(sub_index) + ".fq"
                directory = os.path.dirname(os.path.realpath(__file__))
                p = subprocess.Popen([directory+"/renamer.py", "--bam", args.bam, "--barcodes", args.barcodes, "--out", fq_name,
                        "--chrom", chrom, "--start", str(start), "--end", str(end), "--no_umi", str(args.no_umi), 
                        "--umi_tag", args.umi_tag, "--cell_tag", args.cell_tag])
                all_fastqs.append(fq_name)
                procs[index] = p
                region_fastqs[index].append(fq_name)
                any_running = True
        time.sleep(0.5)
    with open(args.out_dir + "/fastqs.done", 'w') as done:
        for fastqs in region_fastqs:
            done.write("\t".join(fastqs) + "\n")
    return((region_fastqs, all_fastqs))

def remap(args, region_fastqs, all_fastqs):
    print("remapping")
    # run minimap2
    minimap_tmp_files = []
    for index in range(args.threads):
        if index > len(region_fastqs) or len(region_fastqs[index]) == 0:
            continue
        output = args.out_dir + "/souporcell_minimap_tmp_" + str(index) + ".sam"
        minimap_tmp_files.append(output)
        with open(args.out_dir + "/tmp.fq", 'w') as tmpfq:
            subprocess.check_call(['cat'] + region_fastqs[index], stdout = tmpfq)
        with open(output, 'w') as samfile:
            with open(args.out_dir + "/minimap.err",'w') as minierr:
                minierr.write("mapping\n")
                if args.aligner == "HISAT2":
                    fasta_base = args.fasta
                    if fasta_base[-3:] == ".fa":
                        fasta_base = fasta_base[:-3]
                    elif fasta_base[-6:] == ".fasta":
                        fasta_base = fasta_base[:-6]
                    else:
                        assert False, "fasta file not right extension .fa or .fasta"
                    subprocess.check_call(["hisat2", "-p", str(args.threads), "-q", args.out_dir + "/tmp.fq", "-x", 
                    fasta_base,
                    "-S", output], stderr =minierr)
                else:
                    cmd = ["minimap2", "-ax", "splice", "-t", str(args.threads), "-G50k", "-k", "21",
                        "-w", "11", "--sr", "-A2", "-B8", "-O12,32", "-E2,1", "-r200", "-p.5", "-N20", "-f1000,5000",
                        "-n2", "-m20", "-s40", "-g2000", "-2K50m", "--secondary=no", args.fasta, args.out_dir + "/tmp.fq"]
                    minierr.write(" ".join(cmd)+"\n")
                    subprocess.check_call(["minimap2", "-ax", "splice", "-t", str(args.threads), "-G50k", "-k", "21", 
                        "-w", "11", "--sr", "-A2", "-B8", "-O12,32", "-E2,1", "-r200", "-p.5", "-N20", "-f1000,5000",
                        "-n2", "-m20", "-s40", "-g2000", "-2K50m", "--secondary=no", args.fasta, args.out_dir + "/tmp.fq"], 
                        stdout = samfile, stderr = minierr)
        subprocess.check_call(['rm', args.out_dir + "/tmp.fq"])

    with open(args.out_dir + '/remapping.done', 'w') as done:
        for fn in minimap_tmp_files:
            done.write(fn + "\n")
    print("cleaning up tmp fastqs")
    # clean up tmp fastqs
    for fq in all_fastqs:
        subprocess.check_call(["rm", fq])
    return(minimap_tmp_files)

def retag(args, minimap_tmp_files):
    print("repopulating cell barcode and UMI tags")
    # run retagger
    procs = []
    retag_files = []
    retag_error_files = []
    retag_out_files = []
    for index in range(args.threads):
        if index > len(minimap_tmp_files) -1:
            continue
        
        outfile = args.out_dir + "/souporcell_retag_tmp_" + str(index) + ".bam"
        retag_files.append(outfile)
        errfile = open(outfile+".err",'w')
        outfileout = open(outfile+".out",'w')
        retag_error_files.append(errfile)
        retag_out_files.append(outfileout)
        print(args.no_umi)
        print(str(args.no_umi))
        cmd = ["retag.py", "--sam", minimap_tmp_files[index], "--no_umi", str(args.no_umi),
            "--umi_tag", args.umi_tag, "--cell_tag", args.cell_tag, "--out", outfile]
        print(" ".join(cmd))
        directory = os.path.dirname(os.path.realpath(__file__))
        p = subprocess.Popen([directory+"/retag.py", "--sam", minimap_tmp_files[index], "--no_umi", str(args.no_umi), 
            "--umi_tag", args.umi_tag, "--cell_tag", args.cell_tag, "--out", outfile], stdout = outfileout, stderr = errfile)
        procs.append(p)
    for (i, p) in enumerate(procs): # wait for processes to finish
        p.wait()
        retag_error_files[i].close()
        retag_out_files[i].close()
        assert not(p.returncode), "retag subprocess ended abnormally with code " + str(p.returncode)
    for outfile in retag_files:
        subprocess.check_call(['rm',outfile+".err", outfile+".out"])


    print("sorting retagged bam files")
    # sort retagged files
    sort_jobs = []
    filenames = []
    with open(args.out_dir + "/retag.err", 'w') as retagerr:
        for index in range(args.threads):
            if index > len(retag_files) - 1:
                continue
            filename = args.out_dir + "/souporcell_retag_sorted_tmp_" + str(index) + ".bam"
            filenames.append(filename)
            p = subprocess.Popen(["samtools", "sort", retag_files[index], '-o', filename], stderr = retagerr)
            sort_jobs.append(p)
        
    # wait for jobs to finish
    for job in sort_jobs:
        job.wait()
        assert not(job.returncode), "samtools sort ended abnormally with code " + str(job.returncode)

    #clean up unsorted bams
    for bam in retag_files:
        subprocess.check_call(["rm", bam])

    print("merging sorted bams")
    final_bam = args.out_dir + "/souporcell_minimap_tagged_sorted.bam"
    subprocess.check_call(["samtools", "merge", final_bam] + filenames)
    subprocess.check_call(["samtools", "index", final_bam])
    
    print("cleaning up tmp samfiles")
    # clean up tmp samfiles
    for samfile in minimap_tmp_files:
        subprocess.check_call(["rm", samfile])

    # clean up tmp bams
    for filename in filenames:
        subprocess.check_call(['rm', filename])
    subprocess.check_call(["touch", args.out_dir + "/retagging.done"])

def freebayes(args, bam, fasta):
    if not(args.common_variants == None) or not(args.known_genotypes == None):
        if not(args.common_variants == None):
            print("using common variants")
        if not(args.known_genotypes == None):
            print("using known genotypes")
            args.common_variants = args.known_genotypes
        
        # parallelize the samtools depth call. It takes too long
        regions = get_bam_regions(bam, int(args.threads))
        depth_files = []
        depth_procs = []
        print(len(regions))
        for (index, region) in enumerate(regions):
            region_args = []
            for (chrom, start, stop) in region:
                region_args.append(chrom+":"+str(start)+"-"+str(stop))
            depthfile = args.out_dir+"/depth_"+str(index)+".bed"
            depth_files.append(depthfile)
            min_cov = int(args.min_ref)+int(args.min_alt)
            with open(depthfile, 'w') as bed:
                with open(depthfile+".sh",'w') as depther:
                    depther.write("samtools view -hb "+bam+" "+" ".join(region_args)+ " | samtools depth - | "+
                    "awk '{ if ($3 >= "+str(min_cov)+ " && $3 < 100000) { print $1 \"\t\" $2 \"\t\" $2+1 \"\t\" $3 } }'")
                subprocess.check_call(["chmod", "777", depthfile+".sh"])
                #ps0 = subprocess.Popen(['samtools', 'view', bam]+region_args, stdout = subprocess.PIPE)
                #ps1 = subprocess.Popen(['samtools', 'depth', '-'], stdin = ps0.stdout, stdout = subprocess.PIPE)
                # awk magic 
                #ps2 = subprocess.Popen(["awk '{ if ($3 >= " + str(min_cov) + " && $3 < 100000) { print $1 \"\t\" $2 \"\t\" $2+1 \"\t\" $3 } }'"], 
                #    shell = True, stdin = ps1.stdout, stdout = bed)
                ps = subprocess.Popen([depthfile+".sh"], shell = True, stdout = bed)
                depth_procs.append(ps)

        for proc in depth_procs:
            proc.wait()
        merged_depthfiles = []
        for depth_file in depth_files:
            merged_depthfile = depth_file[:-4]+"_merged.bed"
            with open(merged_depthfile, 'w') as bed:
                subprocess.check_call(["bedtools", "merge", "-i", depth_file], stdout = bed)
            merged_depthfiles.append(merged_depthfile)
        with open(args.out_dir + "/depth_merged.bed", 'w') as merged_bed:
            subprocess.check_call(['cat']+merged_depthfiles, stdout = merged_bed)
        for tmp in depth_files: # clean up tmp bed files
            subprocess.check_call(['rm', tmp, tmp+".sh"])
        for tmp in merged_depthfiles:
            subprocess.check_call(['rm', tmp])


        with open(args.out_dir + "/common_variants_covered_tmp.vcf", 'w') as vcf:
            subprocess.check_call(["bedtools", "intersect", "-wa", "-a", args.common_variants, "-b", args.out_dir + "/depth_merged.bed"], stdout = vcf)
        with open(args.out_dir + "/common_variants_covered_tmp.vcf") as vcf:
            with open(args.common_variants) as common:
                with open(args.out_dir + "/common_variants_covered.vcf",'w') as out:
                    for line in common:
                        if line.startswith("#"):
                            out.write(line)
                        else:
                            break
                    for line in vcf:
                        out.write(line)
        with open(args.out_dir + "/variants.done", 'w') as done:
            done.write(args.out_dir + "/common_variants_covered.vcf" + "\n")
        return(args.out_dir + "/common_variants_covered.vcf")

    regions = get_fasta_regions(args.fasta, int(args.threads))
    print(regions)

    region_vcfs = [[] for x in range(args.threads)]
    all_vcfs = []
    bed_files = []
    procs = [None for x in range(args.threads)]
    any_running = True
    filehandles = []
    errhandles = []
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
                    assert not(procs[index].returncode), "freebayes subprocess terminated abnormally with code " + str(procs[index].returncode)
            if len(region_vcfs[index]) == len(region):
                block = True
            if not block:
                sub_index = len(region_vcfs[index])
                chrom = region[sub_index][0]
                start = region[sub_index][1]
                end = region[sub_index][2]
                vcf_name = args.out_dir + "/souporcell_" + str(index) + "_" + str(sub_index) + ".vcf"
                filehandle = open(vcf_name, 'w')
                filehandles.append(filehandle)
                errhandle = open(vcf_name + ".err", 'w')
                errhandles.append(errhandle)
                    
                cmd = ["freebayes", "-f", args.fasta, "-iXu", "-C", "2",
                    "-q", "20", "-n", "3", "-E", "1", "-m", "30", 
                    "--min-coverage", str(int(args.min_alt)+int(args.min_ref)), "--pooled-continuous", "--skip-coverage", "100000"]
                
                cmd.extend(["-r", chrom + ":" + str(start) + "-" + str(end)])
                print(" ".join(cmd))
                cmd.append(bam)
                errhandle.write(" ".join(cmd) + "\n")
                p = subprocess.Popen(cmd, stdout = filehandle, stderr = errhandle)
                all_vcfs.append(vcf_name)
                procs[index] = p
                region_vcfs[index].append(vcf_name)
                any_running = True
        time.sleep(1)
    for filehandle in filehandles:
        filehandle.close()
    for errhandle in errhandles:
        errhandle.close()
    print("merging vcfs")
    subprocess.check_call(["ls "+args.out_dir+'/*.vcf | xargs -n1 -P'+str(args.threads)+' bgzip'],shell=True)
    all_vcfs = [vcf+".gz" for vcf in all_vcfs]
    subprocess.check_call(["ls "+args.out_dir+"/*.vcf.gz | xargs -n1 -P"+str(args.threads) +" bcftools index"],shell=True)
    with open(args.out_dir + "/souporcell_merged_vcf.vcf", 'w') as vcfout:
        subprocess.check_call(["bcftools", "concat", '-a'] + all_vcfs, stdout = vcfout)
    with open(args.out_dir + "/bcftools.err", 'w') as vcferr:
        with open(args.out_dir + "/souporcell_merged_sorted_vcf.vcf", 'w') as vcfout:
            subprocess.check_call(['bcftools', 'sort', args.out_dir + "/souporcell_merged_vcf.vcf"], stdout = vcfout, stderr = vcferr)
    if not args.common_variants == None:
        with open(args.out_dir + "/common.err", 'w') as err:
            with open(args.out_dir + "/vcftmp", 'w') as out:
                subprocess.check_call(['bedtools', 'intersect', '-wa', 
                    '-a', args.out_dir + "/souporcell_merged_vcf.vcf", '-b', args.common_variants], stdout = out, stderr = err)
        subprocess.check_call(['mv', args.out_dir + "/vcftmp", args.out_dir + "/souporcell_merged_sorted_vcf.vcf"])
    subprocess.check_call(['rm', args.out_dir + '/souporcell_merged_vcf.vcf'])
    subprocess.check_call(['bgzip', args.out_dir + "/souporcell_merged_sorted_vcf.vcf"])
    final_vcf = args.out_dir + "/souporcell_merged_sorted_vcf.vcf.gz"
    subprocess.check_call(['tabix', '-p', 'vcf', final_vcf])
    for vcf in all_vcfs:
        subprocess.check_call(['rm', vcf[:-3] + ".err"])       
        subprocess.check_call(['rm', vcf +".csi"])
    subprocess.check_call(['rm'] + all_vcfs)
    if len(bed_files) > 0:
        for bed in bed_files:
            subprocess.check_call(['rm', bed + ".bed"])
        subprocess.check_call(['rm'] + bed_files)
        
    with open(args.out_dir + "/variants.done", 'w') as done:
        done.write(final_vcf + "\n")
    return(final_vcf)


def vartrix(args, final_vcf, final_bam):
    print("running vartrix")
    ref_mtx = args.out_dir + "/ref.mtx"
    alt_mtx = args.out_dir + "/alt.mtx"  
    barcodes = args.barcodes
    if barcodes[-3:] == ".gz":
        with open(args.out_dir + "/barcodes.tsv",'w') as bcsout:
            subprocess.check_call(['gunzip', '-c', barcodes],stdout = bcsout)
        barcodes = args.out_dir + "/barcodes.tsv"
    with open(args.out_dir + "/vartrix.err", 'w') as err:
        with open(args.out_dir + "/vartrix.out", 'w') as out:
            cmd = ["vartrix", "--mapq", "30", "-b", final_bam, "-c", barcodes, "--scoring-method", "coverage", "--threads", str(args.threads),
                "--ref-matrix", ref_mtx, "--out-matrix", alt_mtx, "-v", final_vcf, "--fasta", args.fasta]
            if not(args.no_umi) and args.umi_tag == "UB":
                cmd.append("--umi")
            subprocess.check_call(cmd, stdout = out, stderr = err)
    subprocess.check_call(['touch', args.out_dir + "/vartrix.done"])
    subprocess.check_call(['rm', args.out_dir + "/vartrix.out", args.out_dir + "/vartrix.err"])
    return((ref_mtx, alt_mtx))

def souporcell(args, ref_mtx, alt_mtx, final_vcf):
    print("running souporcell clustering")
    cluster_file = args.out_dir + "/clusters_tmp.tsv"
    with open(cluster_file, 'w') as log:
        with open(args.out_dir+"/clusters.err",'w') as err:
            directory = os.path.dirname(os.path.realpath(__file__))
            cmd = [directory+"/souporcell/target/release/souporcell", "-k",args.clusters, "-a", alt_mtx, "-r", ref_mtx, 
                "--restarts", str(args.restarts), "-b", args.barcodes, "--min_ref", args.min_ref, "--min_alt", args.min_alt, 
                "--threads", str(args.threads)]
            if not(args.known_genotypes == None):
                cmd.extend(['--known_genotypes', final_vcf])
                if not(args.known_genotypes_sample_names == None):
                    cmd.extend(['--known_genotypes_sample_names'] + args.known_genotypes_sample_names)
            print(" ".join(cmd))
            subprocess.check_call(cmd, stdout = log, stderr = err) 
    subprocess.check_call(['touch', args.out_dir + "/clustering.done"])
    return(cluster_file)

def doublets(args, ref_mtx, alt_mtx, cluster_file):
    print("running souporcell doublet detection")
    doublet_file = args.out_dir + "/clusters.tsv"
    with open(doublet_file, 'w') as dub:
        with open(args.out_dir+"/doublets.err",'w') as err:
            directory = os.path.dirname(os.path.realpath(__file__))
            subprocess.check_call([directory+"/troublet/target/release/troublet", "--alts", alt_mtx, "--refs", ref_mtx, "--clusters", cluster_file], stdout = dub, stderr = err)
    subprocess.check_call(['touch', args.out_dir + "/troublet.done"])
    return(doublet_file)

def consensus(args, ref_mtx, alt_mtx, doublet_file):
    print("running co inference of ambient RNA and cluster genotypes")
    directory = os.path.dirname(os.path.realpath(__file__))
    subprocess.check_call([directory+"/consensus.py", "-c", doublet_file, "-a", alt_mtx, "-r", ref_mtx, "-p", args.ploidy,
        "--output_dir",args.out_dir,"--soup_out", args.out_dir + "/ambient_rna.txt", "--vcf_out", args.out_dir + "/cluster_genotypes.vcf", "--vcf", final_vcf])
    subprocess.check_call(['touch', args.out_dir + "/consensus.done"])



#### MAIN RUN SCRIPT
if os.path.isdir(args.out_dir):
    print("restarting pipeline in existing directory " + args.out_dir)
else:
    subprocess.check_call(["mkdir", "-p", args.out_dir])
if not args.skip_remap:
    if not os.path.exists(args.out_dir + "/fastqs.done"):
        (region_fastqs, all_fastqs) = make_fastqs(args)
    else:
        all_fastqs = []
        region_fastqs = []
        with open(args.out_dir + "/fastqs.done") as fastqs:
            for line in fastqs:
                toks = line.strip().split("\t")
                region_fastqs.append(toks)
                for tok in toks:
                    all_fastqs.append(tok)
    if not os.path.exists(args.out_dir + "/remapping.done"):
        minimap_tmp_files = remap(args, region_fastqs, all_fastqs)
    else:
        minimap_tmp_files = []
        with open(args.out_dir + "/remapping.done") as bams:
            for line in bams:
                minimap_tmp_files.append(line.strip())
    if not os.path.exists(args.out_dir + "/retagging.done"):
        retag(args, minimap_tmp_files)
    bam = args.out_dir + "/souporcell_minimap_tagged_sorted.bam" 
else:
    bam = args.bam
if not os.path.exists(args.out_dir + "/variants.done"):
    final_vcf = freebayes(args, bam, fasta)
else:
    with open(args.out_dir + "/variants.done") as done:
        final_vcf = done.readline().strip()
if not os.path.exists(args.out_dir + "/vartrix.done"):
    vartrix(args, final_vcf, bam)
ref_mtx = args.out_dir + "/ref.mtx"
alt_mtx = args.out_dir + "/alt.mtx"
if not(os.path.exists(args.out_dir + "/clustering.done")):
    souporcell(args, ref_mtx, alt_mtx, final_vcf)
cluster_file = args.out_dir + "/clusters_tmp.tsv"
if not(os.path.exists(args.out_dir + "/troublet.done")):
    doublets(args, ref_mtx, alt_mtx, cluster_file)
doublet_file = args.out_dir + "/clusters.tsv"
if not(os.path.exists(args.out_dir + "/consensus.done")):
    consensus(args, ref_mtx, alt_mtx, doublet_file)
print("done")

#### END MAIN RUN SCRIPT        


