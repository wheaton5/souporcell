#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(
    description="\\TODO description")
parser.add_argument("-n", "--num_data_points", required = True, type = int, help = "")
parser.add_argument("-k", "--num_clusters", required = True, type = int, help = "")
parser.add_argument("--num_dimensions", required = True, type = int, help = "")
parser.add_argument("-l", "--lambda", required = True, type = float, help = "", dest="lambda_")
parser.add_argument("-o", "--out_dir", required = True, help = "name of directory to place cellector files")
args = parser.parse_args()

import numpy as np
import math


# low and high might need to be changed if something somewhere else breaks
probabilities = [[np.random.uniform(low=0, high=1) for dimension in range(args.num_dimensions)] for cluster in range(args.num_clusters)]

data_points = []
for dp_idx in range(args.num_data_points):
    cluster = math.floor(args.num_clusters * dp_idx / args.num_data_points)
    data_point = []
    for dimension in range(args.num_dimensions):
        prob = probabilities[cluster][dimension]
        num = np.random.poisson(args.lambda_)
        num_alt = np.random.binomial(num, prob)
        num_ref = num - num_alt
        data_point.append((num_alt, num_ref))
    data_points.append(data_point)

# outputs:

with open(args.out_dir+"/barcodes.tsv", "w") as bar_fid:
    for i in range(len(data_points)):
        bar_fid.write(str(i+1)+"\n")

with open(args.out_dir+"/alt.mtx", "w") as alt_fid, open(args.out_dir+"ref.mtx", "w") as ref_fid:
    alt_fid.write("%%MatrixMarket matrix coordinate real general\n")
    alt_fid.write("% randomly generated\n")
    alt_fid.write("{} {} {}\n".format(arg.num_dimensions, args.num_datapoints, 0))  # num_loci num_cells unread
    ref_fid.write("%%MatrixMarket matrix coordinate real general\n")
    ref_fid.write("% randomly generated\n")
    ref_fid.write("{} {} {}\n".format(arg.num_dimensions, args.num_datapoints, 0))
    for dimension in range(args.num_dimensions):
        for dp_idx in range(len(data_points)):
            data_point = data_points[dp_idx]
            num_alt = data_point[dimension][0]
            num_ref = data_point[dimension][1]
            if num_alt != 0 or num_ref != 0:
                alt_fid.write("{} {} {}\n".format(dimension, dp_idx+1, num_alt))
                ref_fid.write("{} {} {}\n".format(dimension, dp_idx+1, num_ref))

# \\TODO gt.tsv

