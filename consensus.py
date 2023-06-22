#!/usr/bin/env python
import pystan
import argparse
import vcf

parser = argparse.ArgumentParser(description="consensus genotype calling and ambient RNA estimation")
parser.add_argument("-c","--clusters",required=True,help="tsv cluster file from the troublet output")
parser.add_argument("-a","--alt_matrix",required=True,help="alt matrix file")
parser.add_argument("-r","--ref_matrix",required=True,help="ref matrix file")
parser.add_argument("-p","--ploidy",required=False, help="ploidy, must be 1 or 2, defaults to 2")
parser.add_argument("--soup_out",required=True, help="soup output")
parser.add_argument("--vcf_out",required=True, help="vcf output")
parser.add_argument("--output_dir", required=True, help = "output directory")
#parser.add_argument("-d","--doublets",required=True, help="doublet calls")
parser.add_argument("-v","--vcf",required=True,help="vcf file from which alt and ref matrix were created")
args = parser.parse_args()


dirname = args.output_dir


def myopen(fname): return open(fname, 'rb') if fname.endswith('.gz') else open(fname)
cell_genotype_consensus = """
data {
    int<lower=0> cells; // number of cells
    int<lower=0> loci; // number of loci
    int<lower=0> msoup; //number of loci used to estimate soup
    int<lower=0> k; // number of clusters
    int<lower=1,upper=2> ploidy;    
    int cluster_allele_counts[loci, k, 2];
    int cluster_allele_counts_soup[msoup,k,2];
    int cluster_num_cells[k];
    real average_allele_expression[loci, 2];
    real average_allele_expression_soup[msoup,2];
}

transformed data {
    real<upper=0> neg_log_3;
    real<upper=0> fp_prior;
    real<upper=0> tp_prior;
    real<upper=0> logp_base_correct;
    neg_log_3 = -log(ploidy+1);
    tp_prior = log(0.9);
    fp_prior = log(0.1);
    logp_base_correct = log(0.99);
}

parameters {
    real<lower=0.0, upper=1.0> p_soup;
}

model {
    real stuff[ploidy+1];
    real p_hom_ref;
    real p_hom_alt;
    real p_het_ref;
    real truth;
    real p_err;
    real err;
    real hom_ref;
    real hom_alt;
    real het;
    real lse;
    int depth;
    p_soup ~ beta(2,8);
    for (locus in 1:msoup){
        err = 0;
        truth = 0;
        for (cluster in 1:k) {
            depth = cluster_allele_counts_soup[locus][cluster][1] + cluster_allele_counts_soup[locus][cluster][2];
            if (depth == 0) {
                continue;
            }
            p_hom_ref = (1.0 - p_soup) * 1.0 + p_soup * average_allele_expression_soup[locus][1]/(average_allele_expression_soup[locus][1] + average_allele_expression_soup[locus][2]);
            p_hom_alt = (1.0 - p_soup) * 1.0 + p_soup * average_allele_expression_soup[locus][2]/(average_allele_expression_soup[locus][1] + average_allele_expression_soup[locus][2]);
            p_het_ref = (1.0 - p_soup) * 0.5 + p_soup * average_allele_expression_soup[locus][1]/(average_allele_expression_soup[locus][1] + average_allele_expression_soup[locus][2]);
            p_err =  average_allele_expression_soup[locus][1]/(average_allele_expression_soup[locus][1] + average_allele_expression_soup[locus][2]);
            
            hom_ref = binomial_lpmf(cluster_allele_counts_soup[locus][cluster][1] | depth, p_hom_ref);

            hom_alt = binomial_lpmf(cluster_allele_counts_soup[locus][cluster][2] | depth, p_hom_alt);
            
            het = binomial_lpmf(cluster_allele_counts_soup[locus][cluster][1] | depth, p_het_ref);
            
            err += binomial_lpmf(cluster_allele_counts_soup[locus][cluster][1] | depth, p_err);

            stuff[1] = hom_ref + neg_log_3;
            stuff[2] = hom_alt + neg_log_3;
            //stuff2 += err + fp_prior;
            if (ploidy == 2) {
                stuff[3] = het + neg_log_3;
            }
            truth += log_sum_exp(stuff);
        }
        target += log_mix(0.5, truth, err);
    }
}
generated quantities {
    real genotypes[loci,k,ploidy+1];
    real stuff;
    real p_hom_ref;
    real truth[loci];
    real p_hom_alt;
    real p_het_ref;
    real p_err;
    real err[loci];
    real hom_ref;
    real hom_alt;
    real het;
    real lse;
    int depth;
    for (locus in 1:loci){
        err[locus] = 0;
        truth[locus] = 0;
        for (cluster in 1:k) {
            depth = cluster_allele_counts[locus][cluster][1] + cluster_allele_counts[locus][cluster][2];
            if (depth == 0) {
                genotypes[locus][cluster][1] = neg_log_3;
                genotypes[locus][cluster][2] = neg_log_3;
                if (ploidy == 2) {
                    genotypes[locus][cluster][3] = neg_log_3;
                }
                truth[locus] += log_sum_exp(genotypes[locus][cluster]);
                continue;
            }
            p_hom_ref = (1.0 - p_soup) * 1.0 + p_soup * average_allele_expression[locus][1]/(average_allele_expression[locus][1] + average_allele_expression[locus][2]);
            p_hom_alt = (1.0 - p_soup) * 1.0 + p_soup * average_allele_expression[locus][2]/(average_allele_expression[locus][1] + average_allele_expression[locus][2]);
            p_het_ref = (1.0 - p_soup) * 0.5 + p_soup * average_allele_expression[locus][1]/(average_allele_expression[locus][1] + average_allele_expression[locus][2]);
            p_hom_ref = fmax(.01,fmin(.99, p_hom_ref));
            p_hom_alt = fmax(.01,fmin(.99, p_hom_alt));
            p_het_ref = fmax(.01,fmin(.99, p_het_ref));
            p_err =  average_allele_expression[locus][1]/(average_allele_expression[locus][1] + average_allele_expression[locus][2]);
            p_err = fmax(.01,fmin(.99,p_err));
            hom_ref = binomial_lpmf(cluster_allele_counts[locus][cluster][1] | depth, p_hom_ref);
            hom_alt = binomial_lpmf(cluster_allele_counts[locus][cluster][2] | depth, p_hom_alt);
            het = binomial_lpmf(cluster_allele_counts[locus][cluster][1] | depth, p_het_ref);
            err[locus] += neg_log_3 + depth*logp_base_correct + binomial_lpmf(cluster_allele_counts[locus][cluster][1] | depth, p_err);
            genotypes[locus][cluster][1] = hom_ref + neg_log_3;
            genotypes[locus][cluster][2] = hom_alt + neg_log_3;
            if (ploidy == 2) {
                genotypes[locus][cluster][3] = het + neg_log_3;
            }
            truth[locus] += log_sum_exp(genotypes[locus][cluster]);
        }
        truth[locus] += tp_prior;
        err[locus] += fp_prior;
    }
}
"""

if args.ploidy:
    assert(int(args.ploidy) == 1 or int(args.ploidy)==2)
else:
    args.ploidy = 2
import os
import pickle
#sm = pickle.load(open(os.path.realpath(__file__)[0:-12]+"stan_consensus.pickle",'rb'))
sm = pystan.StanModel(model_code=cell_genotype_consensus)
with open("stan_consensus.pickle",'wb') as model:
    pickle.dump(sm, model)
from scipy.special import logsumexp
vcftemplate = vcf.Reader(myopen(args.vcf))
potential_RNAedits = set()
excluded = 0
for (index, rec) in enumerate(vcftemplate):
    if (rec.REF == "T" and str(rec.ALT[0]) == "C") or (rec.REF == "A" and str(rec.ALT[0]) == "G"):
        potential_RNAedits.add(index+1)
        excluded += 1
print(str(excluded) +" excluded for potential RNA editing")


#print("got here")

doublets = set()
with open(args.clusters) as dubs:
    dubs.readline() # get rid of header
    for (index, line) in enumerate(dubs):
        if "doublet" in line or "unassigned" in line:
            doublets.add(index)
print(str(len(doublets))+ " doublets excluded from genotype and ambient RNA estimation")
import subprocess


cell_clusters = {}
cluster_cells = {}
max_cluster = -1
total_cells = 0
cell_size = {}
with open(args.clusters) as assignments:
    assignments.readline()
    for (index, line) in enumerate(assignments):
        if index in doublets:
            continue
        tokens = line.strip().split()
        cell = tokens[0]
        #    doublets.add(index)
        #    continue
        total_cells += 1
        cluster = int(tokens[2])
        cell_clusters[index] = cluster
        cluster_cells.setdefault(cluster,[])
        cluster_cells[cluster].append(index)
        max_cluster = max(max_cluster,cluster)

cluster_counts = [0]*(max_cluster+1)
for (cell, cluster) in cell_clusters.items():
    cluster_counts[cluster] += 1
cell_index = {}
total_lost = 0
loci_counts = {}
loci_full_counts = {}
cell_counts = {}

with open(args.ref_matrix) as ref:
    with open(args.alt_matrix) as alt:
        alt.readline()
        alt.readline()
        tokens = alt.readline().strip().split()
        cells = int(tokens[1])
        total_loci = int(tokens[0])
        ref.readline()
        ref.readline()
        ref.readline()
        for (refline, altline) in zip(ref,alt):
            reftokens = refline.strip().split()
            alttokens = altline.strip().split()
            locus = int(reftokens[0])
            cell = int(reftokens[1])
            if cell - 1 in doublets:
                continue
            cell_counts.setdefault(cell,{})
            cell_counts[cell].setdefault(locus,[0,0])
            loci_counts.setdefault(locus,[0,0])
            refcount = int(reftokens[2])
            altcount = int(alttokens[2])
            cell_counts[cell][locus][0] = refcount
            cell_counts[cell][locus][1] = altcount
            loci_counts.setdefault(locus,[0,0])
            loci_full_counts.setdefault(locus,[0,0])
            loci_full_counts[locus][0] += refcount
            loci_full_counts[locus][1] += altcount
            if refcount > 0:
                loci_counts[locus][0] += 1
            if altcount > 0:
                loci_counts[locus][1] += 1
total_cells = cells
loci_for_soup = {}
index = 0
min_value = 20
excluded = 0
for (locus, counts) in sorted(loci_full_counts.items()):
    if counts[0] > min_value and counts[1] > min_value:
        if not(locus in potential_RNAedits):
            loci_for_soup[locus-1] = index
            index += 1
        else:
            excluded += 1
print(str(excluded)+ " not used for soup calculation due to possible RNA edit")
min_ref = 0
min_alt = 1
used_loci = []
locus_index = {}
index = 0
for (locus, counts) in loci_counts.items():
    if counts[0] >= min_ref and counts[1] >= min_alt:
        used_loci.append(locus-1)
        locus_index[locus-1] = index
        index += 1
used_loci = sorted(used_loci)
used_loci_indices = {locus:i for (i, locus) in enumerate(used_loci)}
total_loci = len(used_loci)
stats_cell_loci = []
stats_cell_counts = []
stats_locus_cells = {}
for (cell, counts) in cell_counts.items():
    stats_cell_loci.append(0)
    stats_cell_counts.append(0)
    for (locus, count) in counts.items():
        if (locus-1) in locus_index:
            stats_locus_cells.setdefault(locus,[0,0])
            if count[0] > 0:
                stats_locus_cells[locus][0] += 1
            if count[1] > 0:
                stats_locus_cells[locus][1] += 1
            stats_cell_loci[-1] += 1
            stats_cell_counts[-1] += (count[0]+count[1])

#with open("cell_loci.csv",'w') as cl:
#    cl.write("loci_per_cell\n")
#    for loci in stats_cell_loci:
#        cl.write(str(loci)+"\n")
#with open("loci_cells.csv",'w') as cl:
#    cl.write("cells_per_locus_ref,cells_per_locus_alt\n")
#    for (locus, cells) in stats_locus_cells.items():
#        cl.write(str(cells[0])+","+str(cells[1])+"\n")
cluster_allele_counts = [[[0,0] for c in range(max_cluster+1)] for x in range(total_loci)]
cluster_allele_counts_soup = [[[0,0] for c in range(max_cluster+1)] for x in range(len(loci_for_soup))]
average_allele_expression_soup = [[0,0] for c in range(len(loci_for_soup))]
average_allele_expression = [[0,0] for c in range(total_loci)]
for (cell, loci_counts) in cell_counts.items():
    if cell-1 in doublets:
        continue
    cluster = cell_clusters[cell-1]
    for (locus, counts) in loci_counts.items():
        ref = 0
        alt = 0
        if counts[0] > 0:
            ref = counts[0]
        if counts[1] > 0:
            alt = counts[1]
        if locus-1 in loci_for_soup:
            cluster_allele_counts_soup[loci_for_soup[locus-1]][cluster][0] += ref
            cluster_allele_counts_soup[loci_for_soup[locus-1]][cluster][1] += alt
            average_allele_expression_soup[loci_for_soup[locus-1]][0] += ref
            average_allele_expression_soup[loci_for_soup[locus-1]][1] += alt
        if not locus-1 in locus_index:
            continue
        cluster_allele_counts[locus_index[locus-1]][cluster][0] += ref
        cluster_allele_counts[locus_index[locus-1]][cluster][1] += alt
        average_allele_expression[locus_index[locus-1]][0] += ref
        average_allele_expression[locus_index[locus-1]][1] += alt

for locus in range(len(average_allele_expression)):
    average_allele_expression[locus][0] /= float(total_cells)
    average_allele_expression[locus][1] /= float(total_cells)

for locus in range(len(average_allele_expression_soup)):
    
    average_allele_expression_soup[locus][0] /= float(total_cells)
    average_allele_expression_soup[locus][1] /= float(total_cells)
    #print(cluster_allele_counts
 
counts_dat = {'cells': total_cells,
              'loci': len(cluster_allele_counts),
              'k': max_cluster + 1,
              'cluster_allele_counts': cluster_allele_counts,
              'cluster_num_cells': cluster_counts,
              'ploidy': int(args.ploidy),
              'msoup': len(cluster_allele_counts_soup),
              'cluster_allele_counts_soup':cluster_allele_counts_soup,
              'average_allele_expression_soup':average_allele_expression_soup,
              'average_allele_expression': average_allele_expression}
#cluster_allele_counts = cluster_allele_counts_soup
#locus_index = loci_for_soup
#print("done loading data")

fit = sm.optimizing(data=counts_dat)

#import pickle
#with open("soup.pickle",'wb') as out:
#    pickle.dump(fit,out)
with open(args.soup_out,'w') as soup:
    soup.write("ambient RNA estimated as "+str(float(fit['p_soup'])*100)+"%")
#with open("cluster_genotypes.tsv",'w') as geno:
#    geno.write("chrom\tpos\tref\talt\t"+"\t".join([i for i in range(max_cluster+1)])+"\n")
import math
import numpy as np
import scipy
from collections import namedtuple
CallData = namedtuple('CallData','GT AO RO T E GO GN')
import gzip
vcftemplate = vcf.Reader(myopen(args.vcf))
vcfreader = vcf.Reader(myopen(args.vcf))
import math
tmp_vcf = dirname+"/tempsouporcell.vcf"
with open(tmp_vcf,'w') as geno:
    vcfwriter = vcf.Writer(geno,vcftemplate)
    samples = [str(cluster) for cluster in range(max_cluster+1)]
    vcfwriter.template.samples = samples
    locus = -1
    for rec in vcfreader:
        locus += 1
        if locus in locus_index:
            newrec = vcf.model._Record(rec.CHROM, rec.POS, rec.ID, rec.REF, rec.ALT, rec.QUAL, rec.FILTER, rec.INFO, 'GT:AO:RO:T:E:GO:GN', {str(x):x for x in range(max_cluster+1)})
            calls = []
            
            for cluster in range(max_cluster+1):
                genotypes = fit['genotypes'][locus_index[locus]][cluster]
                sumexp = logsumexp(genotypes)#[math.exp(x) for x in genotypes])
                #//print(genotypes)
                go = []
                for g in genotypes:
                    if math.isnan(g):
                        go.append('NaN')
                    else:
                        go.append(str(int(g)))
                    
                go = ",".join(go)
                
                #if sumexp == 0:
                #    posteriors = [0.333 for x in range(len(genotypes))]
                #else:
                truth = fit['truth'][locus_index[locus]]
                err = fit['err'][locus_index[locus]]
                total_truth_err = logsumexp([truth,err])
                err = np.exp(err-total_truth_err)
                #print(err)
                logpost = [x - sumexp for x in genotypes]
                posteriors = np.exp(logpost)
                gn = []
                for g in logpost:
                    if math.isnan(g):
                        gn.append("NaN")
                    else:
                        gn.append(str(int(g)))
                gn = ",".join(gn)

                largest = np.argmax(posteriors)
                if err > 0.5:
                    newrec.FILTER = ['BACKGROUND']
                #print(posteriors)
                if len(posteriors) == 3:
                    gt = './.'
                    if posteriors[largest] > 0.5:
                        if largest == 0:
                            gt = '0/0'
                        elif largest == 1:
                            gt = '1/1'
                        elif largest == 2:
                            gt = '0/1'
                elif len(posteriors) == 2:
                    gt = '.'
                    if posteriors[largest] > 0.75:
                        gt = str(largest)

                ao = cluster_allele_counts[locus_index[locus]][cluster][1]
                ro = cluster_allele_counts[locus_index[locus]][cluster][0]
                truth = fit['truth'][locus_index[locus]]
                if not math.isnan(truth):
                    truth = int(truth)
                err = fit['err'][locus_index[locus]]
                if not math.isnan(err):
                    err = int(err)
                calls.append(vcf.model._Call(newrec, str(cluster), CallData(gt, ao, ro, truth, err, go, gn)))
            newrec.samples = calls
            vcfwriter.write_record(newrec)
with open(tmp_vcf) as tmp:
    with open(args.vcf_out,'w') as out:
        for line in tmp:
            if line.startswith("#"):
                if line.startswith("#CHROM\tPOS"):
                    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join([str(clust) for clust in range(max_cluster+1)])+"\n")
                else:
                    out.write(line)
            else:
                out.write(line)
subprocess.check_call(["rm",tmp_vcf])
