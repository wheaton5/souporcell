#!/usr/bin/env python

import argparse
import gzip
import os


parser = argparse.ArgumentParser(description = "determine which clusters are the shared clusters between two souporcell runs with shared samples")
parser.add_argument("-1","--experiment1", required = True, help = "souporcell directory for experiment 1")
parser.add_argument("-2","--experiment2", required = True, help = "souporcell directory for experiment 2")
parser.add_argument("-n","--shared", required= True, type = int, help = "how many samples should these experiments share?")
args = parser.parse_args()

if os.path.isfile(args.experiment1+"/souporcell_merged_sorted_vcf.vcf.gz"):
    vcf1 = gzip.open(args.experiment1+"/souporcell_merged_sorted_vcf.vcf.gz",'rt')
    vcf2 = gzip.open(args.experiment2+"/souporcell_merged_sorted_vcf.vcf.gz",'rt')
elif os.path.isfile(args.experiment1+"/common_variants_covered.vcf"):
    vcf1 = open(args.experiment1+"/common_variants_covered.vcf")
    vcf2 = open(args.experiment2+"/common_variants_covered.vcf")
    

line1 = "#"
while line1.startswith("#"):
    line1 = vcf1.readline()
(chr1, pos1, ref1, alt1) = line1.strip().split()[0:4]
line2 = "#"
while line2.startswith("#"):
    line2 = vcf2.readline()
(chr2, pos2, ref2, alt2) = line2.strip().split()[0:4]
assert chr1 == chr2, "vcfs dont start with same chromosome, are you sure you used the same reference?"
last_chr1 = chr1
last_chr2 = chr2
locus1 = 1
locus2 = 1
locus1_matches = []
locus2_matches = []
locus1_matchset = {}
locus2_matchset = {}
save_loci1 = {locus2: (chr1,pos1)}
save_loci2 = {locus2: (chr2,pos2)}
breakme = False
while not breakme:
    if pos1 == None:
        try:
            (chr1, pos1, ref1, alt1) = vcf1.readline().strip().split()[0:4] 
            locus1 += 1
            save_loci1[locus1] = (chr1, pos1)
        except:
            break
    if pos2 == None:
        try:
            (chr2, pos2, ref2, alt2) = vcf2.readline().strip().split()[0:4]
            locus2 += 1
            save_loci2[locus2] = (chr2, pos2)
        except:
            break
    if chr1 == chr2:
        if pos1 == pos2:
            if alt1 == alt2:
                locus1_matches.append(locus1)
                locus2_matches.append(locus2)
                locus1_matchset[locus1] = len(locus1_matches) - 1
                locus2_matchset[locus2] = len(locus2_matches) - 1
                #print("we have a match chr"+chr1+" "+str(pos1)+" respective loci "+str(locus1)+" "+str(locus2))
            pos1 = None
            pos2 = None
        elif pos1 < pos2:
            pos1 = None
        else:
            pos2 = None
    elif not(chr1 == last_chr1):
        while not(chr2 == chr1):
            try:
                (chr2, pos2, ref2, alt2) = vcf2.readline().strip().split()[0:4]
            except:
                breakme = True
                break
            locus2 += 1
        last_chr1 = chr1
        last_chr2 = chr2
    elif not(chr2 == last_chr2):
        while not(chr1 == chr2):
            try:
                (chr1, pos1, ref1, alt1) = vcf1.readline().strip().split()[0:4]
            except:
                breakme = True
                break
            locus1 += 1
        last_chr1 = chr1
        last_chr2 = chr2
    last_chr1 = chr1
    last_chr2 = chr2

print("locus1 matchset "+str(len(locus1_matchset)))

cell_clusters1 = {}
with open(args.experiment1+"/clusters.tsv") as clust:
    clust.readline()
    for (celldex, line) in enumerate(clust):
        toks = line.strip().split()
        if toks[1] == "singlet":
            cell_clusters1[celldex + 1] = int(toks[2])

print("cell clusters "+str(len(cell_clusters1)))
cell_clusters2 = {}
with open(args.experiment2+"/clusters.tsv") as clust:
    clust.readline()
    for (celldex, line) in enumerate(clust):
        toks = line.strip().split()
        if toks[1] == "singlet":
            cell_clusters2[celldex + 1] = int(toks[2]) 

print("cell clusters "+str(len(cell_clusters2)))
cluster1_locus_counts = {} # cluster to locus to count tuple
total_locus_counts1 = {}
with open(args.experiment1+"/alt.mtx") as alts:
    with open(args.experiment1+"/ref.mtx") as refs:
        alts.readline()
        alts.readline()
        alts.readline()
        refs.readline()
        refs.readline()
        refs.readline()
        for (altline, refline) in zip(alts, refs):
            (locus1, cell1, altcount) = altline.strip().split()
            (locus2, cell2, refcount) = refline.strip().split()
            locus1 = int(locus1)
            locus2 = int(locus2)
            cell1 = int(cell1)
            cell2 = int(cell2)
            assert locus1 == locus2, "ref and alt mtx for experiment 1 dont match?"
            assert cell1 == cell2, "ref and alt mtx for experiment 1 dont match?"
            if locus1 in locus1_matchset:
                if cell1 in cell_clusters1:
                    cluster = cell_clusters1[cell1]
                    matchlocus = locus1_matchset[locus1]
                    if not(cluster in cluster1_locus_counts):
                        cluster1_locus_counts[cluster] = [[0,0] for x in range(len(locus1_matchset))]
                    locus_counts = cluster1_locus_counts[cluster] #cluster1_locus_counts.setdefault(cluster, [[0,0] for x in range(len(locus1_matchset))])
                    counts = locus_counts[matchlocus]
                    counts[0] += int(refcount)
                    counts[1] += int(altcount)
                    total_locus_counts1.setdefault(matchlocus, 0)
                    total_locus_counts1[matchlocus] += int(refcount) + int(altcount)
                    #print(total_locus_counts1[matchlocus])
print("clusters for experiment1 "+str(len(cluster1_locus_counts)))
cluster2_locus_counts = {}
total_locus_counts2 = {}
with open(args.experiment2+"/alt.mtx") as alts:
    with open(args.experiment2+"/ref.mtx") as refs:
        alts.readline()
        alts.readline()
        alts.readline()
        refs.readline()
        refs.readline()
        refs.readline()
        for (altline, refline) in zip(alts, refs):
            (locus1, cell1, altcount) = altline.strip().split()
            (locus2, cell2, refcount) = refline.strip().split()
            locus1 = int(locus1)
            locus2 = int(locus2)
            cell1 = int(cell1)
            cell2 = int(cell2)
            assert locus1 == locus2, "ref and alt mtx for experiment 1 dont match?"
            assert cell1 == cell2, "ref and alt mtx for experiment 1 dont match?"
            if locus1 in locus2_matchset:
                if cell1 in cell_clusters2:
                    cluster = cell_clusters2[cell1]
                    matchlocus = locus2_matchset[locus1]
                    if not(cluster in cluster2_locus_counts):
                        cluster2_locus_counts[cluster] = [[0,0] for x in range(len(locus2_matchset))]
                    locus_counts = cluster2_locus_counts[cluster]#cluster2_locus_counts.setdefault(cluster, [[0,0] for x in range(len(locus2_matchset))])
                    counts = locus_counts[matchlocus]
                    counts[0] += int(refcount)
                    counts[1] += int(altcount)
                    total_locus_counts2.setdefault(matchlocus, 0)
                    #print(matchlocus)
                    total_locus_counts2[matchlocus] += int(refcount) + int(altcount)

print("clusters for experiment2 "+str(len(cluster2_locus_counts)))
            
distances = {}
loci_to_use = set()
for locus in range(len(cluster1_locus_counts[0])):
    has_enough = True
    for cluster in cluster1_locus_counts.keys():
        locus_counts = cluster1_locus_counts[cluster]
        if locus_counts[locus][0] + locus_counts[locus][1] < 8:
            has_enough = False
    for cluster in cluster2_locus_counts.keys():
        locus_counts = cluster2_locus_counts[cluster]
        if locus_counts[locus][0] + locus_counts[locus][1] < 8:
            has_enough = False
    if has_enough:
        loci_to_use.add(locus)

for cluster1 in cluster1_locus_counts.keys():
    locus_counts1 = cluster1_locus_counts[cluster1]
    for cluster2 in cluster2_locus_counts.keys():
        locus_counts2 = cluster2_locus_counts[cluster2]
        loss = 0
        for locus in range(len(locus_counts1)):
            #print("are we")
            locusname = locus
            #if locusname in total_locus_counts1 and total_locus_counts1[locusname] > 20 and locusname in total_locus_counts2 and total_locus_counts2[locusname] > 20:
            if locus in loci_to_use:
                counts1 = locus_counts1[locus]
                counts2 = locus_counts2[locus]
                if counts1[0] + counts1[1] > 0:
                    af1 = counts1[0]/(counts1[0] + counts1[1])
                else:
                    af1 = 1.0
                if counts2[0] + counts2[1] > 0:
                    af2 = counts2[0]/(counts2[0] + counts2[1])
                else:
                    af2 = 1.0
                loss += (af1 - af2)**2
        distances[(cluster1, cluster2)] = loss

distances_sorted = sorted(distances.items(), key=lambda kv: kv[1])
print("experiment1_cluster\texperiment2_cluster\tloss")
for (index, ((cluster1, cluster2), loss)) in enumerate(distances_sorted):
    if index >= args.shared*2:
        break
    print("\t".join([str(cluster1), str(cluster2), str(loss)]))
