cell_genotype_consensus = """
data {
    int<lower=0> n; // number of cells
    int<lower=0> m; // number of loci
    int<lower=0> k; // number of clusters
    
    int cluster_allele_counts[m, k, 2];
    int cluster_num_cells[k];
    real average_allele_expression[m, 2];
}

parameters {
    real<lower=0.0, upper=1.0> theta[m, k, 2];
    real<lower=0.0, upper=1.0> p_soup;
}

model {
    for (cluster in 1:k) {
        for (loci in 1:m){
            real p = 0.0;
            for (allele in 1:2) {
                real lambda_soup = p_soup * average_allele_expression[loci][allele]*cluster_num_cells[cluster];
                real soup = poisson_lpmf(cluster_allele_counts[loci][cluster][allele] | 
                    lambda_soup);
                real lambda_nonsoup = lambda_soup + (1.0-p_soup)* 
                    (average_allele_expression[loci][1] + average_allele_expression[loci][2]) *
                    cluster_num_cells[cluster];
                real nonsoup = poisson_lpmf(cluster_allele_counts[loci][cluster][allele] |
                    lambda_nonsoup);
                target += log_mix(theta[loci][cluster][allele] ,nonsoup, soup);
            }
        }
    }
}
"""

import pystan

sm = pystan.StanModel(model_code=cell_genotype_consensus)
print "done compiling model"

doublets = set()
with open("assignments.tsv") as ass:
    for line in ass:
        tokens = line.strip().split()
        cell = tokens[0]
        c1 = int(tokens[1])
        c2 = int(tokens[2])
        post = float(tokens[4])
        if not c1==c2:
            doublets.add(cell)
        if post < 0.9:
            doublets.add(cell)


cell_clusters = {}
cluster_cells = {}
max_cluster = -1
total_cells = 0
cluster_counts = [0]*6
with open("cell_clusters.csv") as clusters:
    for line in clusters:
        tokens = line.strip().split(",")
        cell = tokens[0]
        if cell in doublets:
            continue
        total_cells += 1
        cluster = int(tokens[1])
        cell_clusters[cell] = cluster
        cluster_cells.setdefault(cluster,[])
        cluster_cells[cluster].append(cell)
        cluster_counts[cluster]+=1
        max_cluster = max(max_cluster,cluster)

cluster_allele_counts = []
average_allele_expression = []
all_alleles = []
with open("calls.tsv") as calls:
    for line in calls:
        tokens = line.strip().split("\t")
        cells = tokens[7].split(";")
        ref_cells = cells[0].split(",")
        alt_cells = cells[1].split(",")
        counts = [[0,0] for x in range(max_cluster+1)]
        chrom = tokens[0]
        pos = tokens[1]
        ref = tokens[3]
        alt = tokens[4]
        all_alleles.append((chrom, pos, ref, alt))
        for ref_cell in ref_cells:
            if ref_cell in doublets:
                continue
            cluster = cell_clusters[ref_cell]
            counts[cluster][0] += 1
        for alt_cell in alt_cells:
            if alt_cell in doublets:
                continue
            cluster = cell_clusters[alt_cell]
            counts[cluster][1] += 1
        cluster_allele_counts.append(counts)
        allele_expression = [float(len(ref_cells))/float(total_cells), float(len(alt_cells))/float(total_cells)]
        average_allele_expression.append(allele_expression)

counts_dat = {'n': total_cells,
              'm': len(cluster_allele_counts),
              'k': max_cluster + 1,
              'cluster_allele_counts': cluster_allele_counts,
              'cluster_num_cells': cluster_counts,
              'average_allele_expression': average_allele_expression}
print "done loading data"

fit = sm.optimizing(data=counts_dat)

import pickle
with open("fit.pickle",'w') as out:
    pickle.dump(fit,out)


print float(fit['p_soup'])
with open("theta.tsv",'w') as theta:
    for index, allele in enumerate(fit['theta']):
        theta.write("\t".join([x for x in all_alleles[index]]+[str(x[0])+","+str(x[1]) for x in fit['theta'][index]]) + "\n")
