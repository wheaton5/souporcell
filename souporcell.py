#!/usr/bin/env python

import numpy as np
import argparse
import tensorflow as tf
import vcf

parser = argparse.ArgumentParser(
    description="single cell RNAseq mixed genotype clustering using sparse mixture model clustering with tensorflow.")
parser.add_argument("-a","--alt_matrix",required=True, help="alt matrix output from vartrix in coverage mode")
parser.add_argument("-r","--ref_matrix",required=True, help="ref matrix output from vartrix in coverage mode")
parser.add_argument("-b","--barcodes",required=True, help="barcodes.tsv file from cellranger")
parser.add_argument("-k","--num_clusters",required=True, help="number of clusters to generate")
parser.add_argument("-l","--max_loci",required=False, help="maximum loci to consider per cell",default=2048)
parser.add_argument("--min_alt",required=False, help="minimum number of cells expressing the alt allele to use the locus for clustering",default=10)
parser.add_argument("--min_ref",required=False, help="minimum number of cells expressing the ref allele to use the locus for clustering",default=10)
parser.add_argument("-t","--threads",required=False, help="number of threads to run on",default=8)
parser.add_argument("-o","--out",required=True,help="output file")
parser.add_argument("--restarts",required=False, default = 15, type = int, help = "number of restarts to clustering")
parser.add_argument("--known_genotypes",required=False, default=None, 
    help="known genotypes population vcf, must be the same vcf used in vartrix (souporcell_pipeline deals with this)")
parser.add_argument("--known_genotypes_sample_names",required=False, default=None, nargs='+', 
    help="known genotypes sample names (list with spaces), sample names must be present in known genotypes vcf")

args = parser.parse_args()

np.random.seed(4) # chosen by dice roll, guaranteed to be random --joke-- https://xkcd.com/221/

min_alt = int(args.min_alt)
min_ref = int(args.min_ref)
K = int(args.num_clusters)

max_loci = int(args.max_loci)

cell_index = {}
total_lost = 0
loci_counts = {}
cell_counts = {}
with open(args.alt_matrix) as alt:
    alt.readline()
    alt.readline()
    tokens = alt.readline().strip().split()
    cells = int(tokens[1])
    cell_counts = {x: {} for x in range(1, cells + 1)}
    total_loci = int(tokens[0])
    for line in alt:
        tokens = line.strip().split()
        locus = int(tokens[0])
        cell = int(tokens[1])
        cell_counts.setdefault(cell, {})
        count = int(tokens[2])
        cell_counts[cell][locus] = [0, count]
        loci_counts.setdefault(locus, [0, 0])
        if count > 0:
            loci_counts[locus][1] += 1
with open(args.ref_matrix) as alt:
    alt.readline()
    alt.readline()
    alt.readline()
    for line in alt:
        tokens = line.strip().split()
        locus = int(tokens[0])
        cell = int(tokens[1])
        count = int(tokens[2])
        cell_counts[cell][locus][0] = count
        loci_counts.setdefault(locus, [0, 0])
        if count > 0:
            loci_counts[locus][0] += 1

used_loci_set = set()
used_loci = []
for (locus, counts) in loci_counts.items():
    if counts[0] >= min_ref and counts[1] >= min_alt:
        used_loci.append(locus - 1)
        used_loci_set.add(locus - 1)
used_loci = sorted(used_loci)
used_loci_indices = {locus:i for (i, locus) in enumerate(used_loci)}
loci = len(used_loci)

if not(args.known_genotypes == None):
    reader = vcf.Reader(open(args.known_genotypes))
    if args.known_genotypes_sample_names == None:
        args.known_genotypes_sample_names = reader.samples
    assert int(args.num_clusters) == len(args.known_genotypes_sample_names), "clusters must equal samples for known_genotypes"
    sample_index = {samp: index for (index, samp) in enumerate(args.known_genotypes_sample_names)}
    sample_genotypes = np.random.random((len(args.known_genotypes_sample_names), loci))#[[np.random.random() for x in range(loci)] for y in range(len(args.known_genotypes_sample_names))]
    filled = 0
    unfilled = 0
    for (index, rec) in enumerate(reader):
        if index in used_loci_set:
            for sample in args.known_genotypes_sample_names:
                gt = rec.genotype(sample).data.GT
                if gt == '0|0' or gt == '0/0':
                    freq = 1.0
                    filled += 1
                elif gt == '0|1' or gt == '1|0' or gt == '0/1':
                    freq = 0.5
                    filled += 1
                elif gt == '1|1' or gt == '1/1':
                    freq = 0.0
                    filled += 1
                else:
                    freq = np.random.random()
                    unfilled += 1
                sampledex = sample_index[sample]

                sample_genotypes[sampledex][used_loci_indices[index]] = freq


print("loci being used based on min_alt, min_ref, and max_loci "+str(loci))
cells = len(cell_counts)
#cells = 0
#for cell in cell_counts.keys():
#    cells = max(cell, cells)

cell_data = np.zeros((cells, max_loci))
cell_loci = np.zeros((cells, max_loci))
weights = np.zeros((cells, max_loci))
for cell in cell_counts.keys():
    index = 0
    single_cell_counts = cell_counts[cell]
    for locus in cell_counts[cell].keys():
        locus_counts = single_cell_counts[locus]
        if loci_counts[locus][0] >= min_ref and loci_counts[locus][1] >= min_alt:
            if index < max_loci:
                ref_c = locus_counts[0]
                alt_c = locus_counts[1]
                if ref_c + alt_c > 0:
                    cell_data[cell - 1][index] = float(ref_c)/float(ref_c + alt_c) if ref_c + alt_c > 0 else 0.0
                    cell_loci[cell - 1][index] = used_loci_indices[locus - 1]
                    weights[cell - 1][index] = 1.0
                    index += 1
                total_lost += 1


data = cell_data
data_loci = cell_loci

rng = np.random
if args.known_genotypes == None:
    phi = tf.get_variable(name = "phi",shape=(loci, K), initializer = tf.initializers.random_uniform(minval = 0, maxval = 1), dtype = tf.float64)
else:
    args.restarts = 1
    init = tf.constant(sample_genotypes.T)
    phi = tf.get_variable(name="phi", initializer = init, dtype = tf.float64)
input_data = tf.placeholder("float64", (cells, max_loci)) #tf.constant("input",np.asmatrix(data))
input_loci = tf.placeholder("int32", (cells, max_loci))
weight_data = tf.placeholder("float64", (cells, max_loci)) #tf.constant("weights",np.asmatrix(weights))
loci_per_cell = tf.placeholder("float64", (cells))
trans = tf.transpose(input_data)
broad_trans = tf.broadcast_to(trans,[K,max_loci,cells])
untrans = tf.transpose(broad_trans)
xtest = untrans-tf.gather(phi,input_loci)
weightshape = tf.transpose(tf.broadcast_to(tf.transpose(weight_data),[K,max_loci,cells]))
weighted = weightshape*xtest
powtest = -tf.pow(weighted,2)
sumtest = tf.reduce_sum(powtest,axis=1)
logsum = tf.reduce_logsumexp(sumtest,axis=1)
cost = -tf.reduce_sum(logsum)
#for i in range(K):
#    for j in range(i+1,K):
#        cost = cost - tf.math.minimum(tf.constant(loci, dtype=tf.float64),tf.reduce_sum(tf.pow(phi[:,i]-phi[:,j],2)))
optimizer = tf.train.AdamOptimizer(learning_rate=0.1).minimize(cost)

post = sumtest#tf.transpose(tf.transpose(sumtest) - tf.reduce_logsumexp(sumtest,axis=1))
repeats = args.restarts
posteriors = []
min_cost = None

threads = int(args.threads)
session_conf = tf.ConfigProto(
      intra_op_parallelism_threads=threads,
      inter_op_parallelism_threads=threads)

for repeat in range(repeats):
    init = tf.global_variables_initializer()
    print("repeat "+str(repeat))
    training_epochs = 1000
    last_cost = None
    with tf.Session(config = session_conf) as sess:
        sess.run(init)
        for epoch in range(training_epochs):
            #if epoch < 70:
            sess.run(optimizer, feed_dict={input_data:data, weight_data:weights, input_loci:data_loci})
            #else:
            #    sess.run(optimizer2, feed_dict={input_data:ref_data, weight_data:weights})
            if epoch % 10 == 0:
                c = sess.run(cost, feed_dict={input_data:data, weight_data:weights, input_loci:data_loci})
                
                print("epoch "+str(epoch)+" "+str(c))
                #if last_cost and ((last_cost-c)/c) < 0.0001:
                if min_cost and last_cost and c > min_cost and (last_cost - c)/(c - min_cost) < 0.005:
                    print("bailing out, too little progress toward minimum so far")
                    break
                if last_cost and last_cost - c < 1:
                    last_cost = None
                    break
                last_cost = c
        if min_cost:
            min_cost = min(min_cost, c)
        else:
            min_cost = c

        #boop = sess.run(sumtest, feed_dict={input_data:data, weight_data:weights, input_loci:data_loci})
        posterior = sess.run(post, feed_dict={input_data:data, weight_data:weights, input_loci:data_loci})
        posteriors.append((c,posterior))


posterior = sorted(posteriors, key=lambda x: x[0])
print(posterior[0])
posterior = posterior[0][1]


print(np.argmax(posterior,axis=1))
cluster_counts = np.zeros(K)
for i in range(len(cluster_counts)):
    print(str(i)+"\t"+str(cluster_counts[i]))
clusters = np.argmax(posterior,axis=1)
cluster_posteriors = posterior[clusters]

barcodes = []
with open(args.barcodes) as bcs:
    for line in bcs:
        barcodes.append(line.strip().split()[0])

with open(args.out,'w') as out:
    for c in range(cells):
        out.write(barcodes[c]+"\t"+str(clusters[c])+"\t"+"\t".join([str(x) for x in posterior[c]]))
        out.write("\n")



