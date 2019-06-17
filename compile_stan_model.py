#!/usr/bin/env python
import pystan


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

import os
import pickle
sm = pystan.StanModel(model_code=cell_genotype_consensus)
with open(os.path.realpath(__file__)[0:-21]+"stan_consensus.pickle",'wb') as model:
    pickle.dump(sm, model)
