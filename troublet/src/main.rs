extern crate fnv;
#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate statrs;
use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;

use hashbrown::{HashMap,HashSet};
use fnv::FnvHasher;
use std::hash::BuildHasherDefault;
type FnvHashMap<T,V> = HashMap<T,V, BuildHasherDefault<FnvHasher>>;
type FnvHashSet<T> = HashSet<T, BuildHasherDefault<FnvHasher>>;

use clap::{App};

use std::cmp::max;
use std::cmp::min;
use statrs::distribution::{Binomial, Discrete, Beta, Continuous};
use statrs::function::{beta,factorial};
use std::f64;

fn main() {
    let params = load_params();
    let (cell_clusters, cluster_losses, num_clusters) = load_clusters(&params.clusters);
    let (cluster_allele_counts, cluster_allele_fractions, soup_allele_fractions, cell_allele_counts, locus_counts) = load_allele_data(&params, &cell_clusters, num_clusters);
    call_doublets(&params, cluster_allele_counts, cell_clusters, cluster_allele_fractions, soup_allele_fractions, cell_allele_counts, cluster_losses, locus_counts); 

}

fn call_doublets(params: &Params, mut cluster_allele_counts: FnvHashMap<(usize,usize),Vec<(f64,f64)>>, 
        cell_clusters: Vec<(usize, usize)>, mut cluster_allele_fractions: FnvHashMap<(usize, usize), Vec<f64>>, 
        soup_allele_fractions: Vec<f64>, cell_allele_counts: Vec<Vec<(usize, u64, u64)>>, 
        cluster_losses: Vec<(String, Vec<f64>)>, total_locus_counts: HashMap<usize, usize>) { // good god i should use more structs
    let num_clusters = cluster_losses[0].1.len();
    print!("barcode\tstatus\tassignment\tsinglet_posterior\tdoublet_posterior\tlog_prob_singleton\tlog_prob_doublet\t");
    for cluster in 0..num_clusters {
        print!("cluster{}",cluster);
        if cluster == num_clusters - 1 {println!();}
        else {print!("\t");}
    }
    let mut all_assignments: Vec<Vec<String>> = Vec::new();
    let mut all_removed: HashSet<usize> = HashSet::new();
    let mut any_removed = 1; // fake for while loop
    let mut best_singlets: Vec<usize> = Vec::new();
    let mut all_cluster_logprobs: Vec<Vec<f64>> = Vec::new();
    while any_removed > 0 {
        all_cluster_logprobs = Vec::new();
        all_assignments = Vec::new();
        best_singlets = Vec::new();
        let mut new_removed: Vec<usize> = Vec::new();
        for (cell, loci) in cell_allele_counts.iter().enumerate() {
            let debug = params.debug.contains(&cell);
            let mut best_doublet1 = 0;
            let mut best_doublet2 = 0;
            let mut best_doublet_log_prob = f64::NEG_INFINITY;
            let mut best_singlet = 0;
            let mut best_singleton_log_prob = f64::NEG_INFINITY;
            let mut singletons: Vec<(f64, usize)> = Vec::new();
            let mut cluster_logprobs: Vec<f64> = Vec::new();
            for cluster1 in 0..num_clusters {
                let singleton_allele_fractions = cluster_allele_fractions.get(&(cluster1, cluster1)).unwrap();
                let mut log_prob_singleton = (1.0 - params.doublet_prior).ln();
                let cluster_counts_vec = cluster_allele_counts.get(&(cluster1,cluster1)).unwrap();
                for (locus, ref_count, alt_count) in loci {
                    if *total_locus_counts.get(locus).unwrap() < 20 { continue; }
                    //if ref_count + alt_count < 2 { continue; }
                    //if !legal_loci[cluster1].contains(&locus) { continue; }
                    let locus_counts = cluster_counts_vec[*locus];
                    //if locus_counts.0 + locus_counts.1 <= 35.0 { continue; }
                    let singleton_allele_fraction = (1.0-params.p_soup)*singleton_allele_fractions[*locus] + params.p_soup*soup_allele_fractions[*locus];
                    let count = locus_counts.0 + locus_counts.1;
                    //if count < 35.0 { continue; }
                    let lbetapdfsingleton = factorial::ln_binomial(*ref_count + *alt_count, *ref_count) + 
                        beta::ln_beta(*ref_count as f64 + singleton_allele_fraction*count + 1.0, 
                            *alt_count as f64 + (1.0 - singleton_allele_fraction)*count + 1.0) - 
                        beta::ln_beta(singleton_allele_fraction*count + 1.0, 
                            (1.0-singleton_allele_fraction)*count + 1.0); //beta.pdf(singleton_allele_fraction).log10();
                    let singleton_logprob = lbetapdfsingleton; //singleton_binomial.pmf(*ref_count).ln();// + lbetapdfsingleton;
                    //let singleton_binomial = Binomial::new(singleton_allele_fraction, ref_count+alt_count).unwrap();
                    log_prob_singleton += singleton_logprob;
                }
                singletons.push((log_prob_singleton, cluster1));
                if log_prob_singleton > best_singleton_log_prob {
                    best_singlet = cluster1;
                    best_singleton_log_prob = log_prob_singleton;
                }
            }
            for (logprob, _) in singletons.iter() {
                cluster_logprobs.push(*logprob);
            }
            singletons.sort_by(|(a, _),(b, _)| (-a).partial_cmp(&-b).unwrap());
            for cluster1 in 0..(num_clusters.min(4)) {
                for cluster2 in (cluster1 + 1)..(num_clusters.min(4)) {
                    let cluster1 = singletons[cluster1].1;
                    let cluster2 = singletons[cluster2].1;

                    let doublet_allele_fractions = cluster_allele_fractions.get(&(cluster1,cluster2)).unwrap();
                    
                    let mut log_prob_doublet = params.doublet_prior.ln();
                    let cluster_counts_vec = cluster_allele_counts.get(&(cluster1,cluster1)).unwrap();
                    let cluster_counts_vec2 = cluster_allele_counts.get(&(cluster1,cluster2)).unwrap();
                    let cluster_counts_vec3 = cluster_allele_counts.get(&(cluster2,cluster2)).unwrap();
                    for (locus, ref_count, alt_count) in loci {

                        if *total_locus_counts.get(locus).unwrap() < 20 { continue; }
                        //if ref_count + alt_count < 2 { continue; }
                        //if !legal_loci[cluster1].contains(&locus) { continue; }
                        let locus_counts = cluster_counts_vec[*locus];
                        //if locus_counts.0 + locus_counts.1 <= 35.0 { continue; }
                        //let beta = Beta::new(1.0+locus_counts.0, 1.0+locus_counts.1).unwrap();
                        let locus_counts2 = cluster_counts_vec2[*locus];
                        let locus_counts3 = cluster_counts_vec3[*locus];
                        
                        //let beta2 = Beta::new(1.0+locus_counts2.0, 1.0+locus_counts2.1).unwrap();
                        
                        let doublet_allele_fraction = (1.0-params.p_soup)*doublet_allele_fractions[*locus] + 
                            params.p_soup*soup_allele_fractions[*locus];
                        //let count = ((locus_counts.0+locus_counts.1)+(locus_counts3.0 + locus_counts3.1))/2.0;
                        let count = ((locus_counts.0+locus_counts.1)+(locus_counts3.0 + locus_counts3.1))/2.0;
                        //if count < 35.0 { continue; }
                        let lbetapdfdoublet = factorial::ln_binomial(*ref_count + *alt_count, *ref_count) + 
                            beta::ln_beta(*ref_count as f64 + doublet_allele_fraction*count + 1.0, 
                                *alt_count as f64 + (1.0-doublet_allele_fraction)*count + 1.0) - 
                            beta::ln_beta(doublet_allele_fraction*count + 1.0, 
                                (1.0-doublet_allele_fraction)*count + 1.0); //beta.pdf(singleton_allele_fraction).log10();
                        //if (singleton_allele_fraction-doublet_allele_fraction).abs() < 0.2 { continue; }
                        
                        //let doublet_binomial = Binomial::new(doublet_allele_fraction, ref_count+alt_count).unwrap();
                        
                        let doublet_logprob = lbetapdfdoublet;//doublet_binomial.pmf(*ref_count).ln();// + lbetapdfdoublet;
                        
                        log_prob_doublet += doublet_logprob;
                    }
                    if log_prob_doublet > best_doublet_log_prob {
                        best_doublet1 = cluster1;
                        best_doublet2 = cluster2;
                        best_doublet_log_prob = log_prob_doublet;
                    }

                }
            }
            let cell_barcode = &cluster_losses[cell].0;
            //let losses: Vec<f64> = cluster_losses[cell].1;
            let mut doublet_question: Vec<f64> = Vec::new();
            //doublet_question.push(log_prob_singleton);
            //doublet_question.push(log_prob_doublet);
            doublet_question.push(best_singleton_log_prob);
            doublet_question.push(best_doublet_log_prob);
            let doublet_question_denom = log_sum_exp(&doublet_question);
            let mut singlet_question: Vec<f64> = Vec::new();
            let mut top_singlet: f64 = f64::NEG_INFINITY;
            let mut singlet_assignment = 0;
            for (cluster, loss) in cluster_losses[cell].1.iter().enumerate() {
                singlet_question.push(*loss);
                if loss > &top_singlet { 
                    top_singlet = *loss; 
                    singlet_assignment = cluster;
                }
            }
            //assert!(singlet_assignment == best_singlet, "{}\t{}\t{}", cell_barcode, singlet_assignment, best_singlet);
            
            best_singlets.push(singlet_assignment);
            let singlet_question_denom = log_sum_exp(&singlet_question);
            
            let singlet_posterior = (top_singlet - singlet_question_denom).exp();
            let doublet_posterior = (best_doublet_log_prob - doublet_question_denom).exp();
            let mut assignment: Vec<String> = Vec::new(); assignment.push(cell_barcode.to_string());
            if best_singleton_log_prob > best_doublet_log_prob {
                if singlet_posterior > params.singlet_threshold {
                    assignment.push(format!("singlet\t{}\t{}\t{}\t{}\t{}", best_singlet, singlet_posterior, doublet_posterior, best_singleton_log_prob, best_doublet_log_prob));
                    //print!("{}\tsinglet\t{}\t{}\t{}\t", cell_barcode, best_singlet, best_singleton_log_prob, best_doublet_log_prob);
                } else {
                    //if !all_removed.contains(&cell) {
                    //    all_removed.insert(cell);
                    //    new_removed.push(cell);
                    //    eprintln!("removing {} as unassigned", cell_barcode);
                    //}
                    assignment.push(format!("unassigned\t{}\t{}\t{}\t{}\t{}", best_singlet, singlet_posterior, doublet_posterior, best_singleton_log_prob, best_doublet_log_prob)); 
                    //print!("{}\tunassigned\t{}\t{}\t{}\t", cell_barcode, best_singlet, best_singleton_log_prob, best_doublet_log_prob);
                }
            } else {
                if doublet_posterior >= params.doublet_threshold {
                    if !all_removed.contains(&cell) {
                        all_removed.insert(cell);
                        new_removed.push(cell);
                        eprintln!("removing {} as doublet", cell_barcode);
                    }
                    
                    assignment.push(format!("doublet\t{}/{}\t{}\t{}\t{}\t{}", best_doublet1, best_doublet2, singlet_posterior, doublet_posterior, best_singleton_log_prob, best_doublet_log_prob));
                    //print!("{}\tdoublet\t{}/{}\t{}\t{}\t", cell_barcode, best_doublet1, best_doublet2, best_singleton_log_prob, best_doublet_log_prob);
                } else {
                    //if !all_removed.contains(&cell) {
                    //    all_removed.insert(cell);
                    //    new_removed.push(cell);
                    //    eprintln!("removing {} as unassigned doublet", cell_barcode);
                   // }
                    assignment.push(format!("unassigned\t{}/{}\t{}\t{}\t{}\t{}", best_doublet1, best_doublet2, singlet_posterior, doublet_posterior, best_singleton_log_prob, best_doublet_log_prob));
                    //print!("{}\tunassigned\t{}/{}\t{}\t{}\t", cell_barcode, best_doublet1, best_doublet2, best_singleton_log_prob, best_doublet_log_prob); 
                }
            }
            all_assignments.push(assignment);
            if singlet_assignment != best_singlet {
                if !all_removed.contains(&cell){
                    eprintln!("error\t{}\t{}\t{}", cell_barcode, singlet_assignment, best_singlet);
                }
            }
            all_cluster_logprobs.push(cluster_logprobs);
            //for (index, loss) in cluster_losses[cell].1.iter().enumerate() {
            //    print!("{}",loss);
            //    if index == num_clusters - 1 {println!();}
            //    else { print!("\t"); }
            //}
        }
        eprintln!("{} cells removed this round", new_removed.len());
        any_removed = new_removed.len();
        for cell in new_removed {
            let cluster1 = best_singlets[cell];
            for cluster2 in 0..num_clusters {
                let mut counts = cluster_allele_counts.get_mut(&(cluster1, cluster2)).unwrap();
                let loci = &cell_allele_counts[cell];
                for (locus, ref_count, alt_count) in loci {
                    counts[*locus].0 = (counts[*locus].0 - *ref_count as f64).max(0.0);
                    counts[*locus].1 = (counts[*locus].1 - *alt_count as f64).max(0.0);
                }
            }
        }
        
        for cluster1 in 0..num_clusters {
            for cluster2 in 0..num_clusters {
            //let cluster2 = cluster1;
                let mut fractions = cluster_allele_fractions.get_mut(&(cluster1,cluster2)).unwrap();
                let cluster_counts = cluster_allele_counts.get_mut(&(cluster1,cluster2)).unwrap();
                for locus in 0..fractions.len() {
                    let counts = cluster_counts[locus];
                    if counts.0 + counts.1 > 0.0 {
                        fractions[locus] = (counts.0/(counts.0+counts.1)).min(0.99).max(0.01);
                    } //else { //eprintln!("0 counts, does this fuck anything up {} ",fractions[locus]); };
                }
            }
        }
    }

    for (cell, assignment) in all_assignments.iter().enumerate() {
        let logprobs = &all_cluster_logprobs[cell];
        for s in assignment {
            print!("{}\t",s);
        }
        for (i, s) in logprobs.iter().enumerate() {
            print!("{}",s);
            if i < logprobs.len() - 1 { print!("\t"); }
        }println!();
    }
}

fn log_sum_exp(p: &Vec<f64>) -> f64{
    let max_p: f64 = p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let sum_rst: f64 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn load_allele_data(params: &Params, cell_clusters: &Vec<(usize, usize)>, num_clusters: usize) -> 
        (FnvHashMap<(usize,usize),Vec<(f64, f64)>>, 
        FnvHashMap<(usize, usize), Vec<f64>>, Vec<f64>, 
        Vec<Vec<(usize, u64, u64)>>, HashMap<usize, usize>) {
    let mut locus_counts: HashMap<usize, usize> = HashMap::new();
    let mut cluster_allele_counts: FnvHashMap<(usize, usize), Vec<(f64, f64)>> = FnvHashMap::default(); 
        // this will be (clust1, clust2) to Vec<f64> which is locus_index to ref allele fraction
    let mut cell_allele_counts: Vec<Vec<(usize, u64, u64)>> = Vec::new(); 
    let mut soup_allele_counts: Vec<(f64, f64)> = Vec::new();
        // so this is going to be cell_index to (locus_index, ref_count, alt_count)
    let alts = File::open(&params.alts).expect("unable to open alt file");
    let alts = BufReader::new(alts);
    let refs = File::open(&params.refs).expect("unable to open ref file");
    let refs = BufReader::new(refs);
    let mut legal_loci: Vec<FnvHashSet<usize>> = Vec::new();
    for cluster in 0..num_clusters {
        let mut lociset: FnvHashSet<usize> = FnvHashSet::default();
        legal_loci.push(lociset);
    }
    for clust1 in 0..num_clusters {
        for clust2 in 0..num_clusters {
            let mut cluster_pair_vector: Vec<(f64, f64)> = Vec::new();
            cluster_allele_counts.insert((clust1, clust2), cluster_pair_vector);
        }
    }
    let mut loci = 0;
    for (index, (refline, altline)) in refs.lines().zip(alts.lines()).enumerate() {
        let refline = refline.expect("could not read line");
        let altline = altline.expect("count not read line");
        if refline.starts_with("%") { continue; }
        if index == 2 {
            let tokens: Vec<&str> = refline.trim().split_whitespace().collect();
            loci = tokens[0].parse::<usize>().unwrap();
            let cells = tokens[1].parse::<usize>().unwrap();
            for clust1 in 0..num_clusters {
                for clust2 in 0..num_clusters {
                    let mut loci_vec = cluster_allele_counts.get_mut(&(clust1, clust2)).unwrap();
                    for _locus in 0..loci {
                        loci_vec.push((0.0, 0.0));
                    }
                }
            }
            for locus in 0..loci {
                soup_allele_counts.push((0.0,0.0));
                locus_counts.insert(locus,0);
            }
            for _cell in 0..cells {
                let mut cell_vec: Vec<(usize, u64, u64)> = Vec::new();
                cell_allele_counts.push(cell_vec);
            }
            continue;
        }
        let ref_tokens: Vec<&str> = refline.trim().split_whitespace().collect();
        let alt_tokens: Vec<&str> = altline.trim().split_whitespace().collect();
        let locus = ref_tokens[0].parse::<usize>().unwrap() - 1;
        let altlocus = alt_tokens[0].parse::<usize>().unwrap() - 1;
        assert!(locus == altlocus,"allele matrices do not match");
        let cell = ref_tokens[1].parse::<usize>().unwrap() - 1;
        let altcell =  alt_tokens[1].parse::<usize>().unwrap() - 1;
        assert!(cell == altcell,"allele matrices do not match");
        let refcount = ref_tokens[2].parse::<u64>().unwrap();
        let altcount = alt_tokens[2].parse::<u64>().unwrap();
        let count = locus_counts.entry(locus).or_insert(0);
        *count += 1;
        let clust1 = cell_clusters[cell].0;
        //let clust2 = cell_clusters[cell].1;
        soup_allele_counts[locus].0 += refcount as f64;
        soup_allele_counts[locus].1 += altcount as f64;
        for clust2 in 0..num_clusters {
            
            let mut locus_vec = cluster_allele_counts.get_mut(&(clust1, clust2)).unwrap();
            locus_vec[locus].0 += refcount as f64;
            locus_vec[locus].1 += altcount as f64;
            if clust1 != clust2 {
                let mut locus_vec = cluster_allele_counts.get_mut(&(clust2, clust1)).unwrap();
                locus_vec[locus].0 += refcount as f64;
                locus_vec[locus].1 += altcount as f64;
            }
        }
        //{
        //    let mut locus_vec2 = cluster_allele_counts.get_mut(&(clust1, clust2)).unwrap();
        //    locus_vec2[locus].0 += refcount as f64;
        //    locus_vec2[locus].1 += altcount as f64;
        //}
        //{   // this scope isn't necessary but included for symmetry
        //    let mut locus_vec3 = cluster_allele_counts.get_mut(&(clust2, clust1)).unwrap();
        //    locus_vec3[locus].0 += refcount as f64;
        //    locus_vec3[locus].1 += altcount as f64;
        //}
        //}
        cell_allele_counts[cell].push((locus, refcount, altcount));
    }
    let mut cluster_allele_fractions: FnvHashMap<(usize, usize), Vec<f64>> = FnvHashMap::default();
    let mut soup_allele_fractions: Vec<f64> = Vec::new();
    for (refcount, altcount) in soup_allele_counts {
        if refcount + altcount > 0.0 {
            soup_allele_fractions.push(refcount/(refcount+altcount));
        } else {
            soup_allele_fractions.push(0.0);
        }
    }
    let min = 0.01;
    let max = 0.99;
    let mut num_0_counts = 0;
    for clust1 in 0..num_clusters {
        for clust2 in 0..num_clusters {
            let mut new_vec: Vec<f64> = Vec::new();
            let counts1 = cluster_allele_counts.get(&(clust1,clust1)).unwrap();
            let counts2 = cluster_allele_counts.get(&(clust2,clust2)).unwrap();
            for locus in 0..loci {
                let counts1 = counts1[locus];
                let counts2 = counts2[locus];
                if counts1.0 + counts1.1 > 35.0 {
                    legal_loci[clust1].insert(locus);
                }
                if counts1.0 + counts1.1 > 0.0 {
                    if counts2.0 + counts2.1 > 0.0 {
                        let mut frac = 0.5*(counts1.0/(counts1.0+counts1.1)) + 0.5*(counts2.0/(counts2.0+counts2.1));
                        if frac > max { frac = max; } else if frac < min { frac = min; }
                        new_vec.push(frac);
                    } else {
                        let mut frac = counts1.0/(counts1.0+counts1.1);
                        if frac > max { frac = max; } else if frac < min { frac = min; }
                        new_vec.push(frac);
                    }
                //} if counts2.0 + counts2.1 > 0.0 {
                //    let mut frac = counts2.0/(counts2.0 + counts2.1);
                //    if frac > max { frac = max; } else if frac < min { frac = min; }
                //    new_vec.push(frac);
                } else {
                    //eprintln!("loaded 0 counts, does this fuck anything up?");
                    num_0_counts += 1;
                    new_vec.push(0.5);
                }
            }
            cluster_allele_fractions.insert((clust1, clust2), new_vec);
        }
    }
    eprintln!("{} loaded 0 counts, is this a problem?", num_0_counts);
    //for (key, loci_vec) in cluster_allele_counts {
    //    let mut new_vec: Vec<f64> = Vec::new();
    //    for (count1, count2) in loci_vec {
    //        if count1 + count2 > 0.0 {
    //            new_vec.push(count1/(count1+count2));
    //        } else {
    //            new_vec.push(0.0);
    //        }
    //    }
    //    cluster_allele_fractions.insert(key, new_vec);
   // }
    (cluster_allele_counts, cluster_allele_fractions, soup_allele_fractions, cell_allele_counts, locus_counts)
}

fn load_clusters(clusters: &String) -> (Vec<(usize, usize)>, Vec<(String, Vec<f64>)>, usize)  {
    let f = File::open(clusters).expect("Unable to open file");
    let f = BufReader::new(f);
    let mut cell_clusters: Vec<(usize, usize)> = Vec::new();
    let mut clusters: Vec<(String, Vec<f64>)> = Vec::new();
    let mut num_clusters: usize = 0;
    for line in f.lines() {
        let line = line.expect("Unable to read line");
        let tokens: Vec<&str> = line.trim().split("\t").collect();
        let cell = tokens[0];
        let mut losses: Vec<f64> = Vec::new();
        let cluster1 = tokens[1].to_string().parse::<usize>().unwrap();
        num_clusters = max(num_clusters, cluster1);
        let mut best: f64 = -10000000000.0;
        let mut second_best: f64 = -10000000000.0;
        let mut second_best_index: usize = 0;
        let mut best_index: usize = 0;
        for i in 2..tokens.len() {
            let value = tokens[i].to_string().parse::<f64>().unwrap();
            losses.push(value);
            if value > best {
                second_best = best;
                second_best_index = best_index;
                best = value;
                best_index = i-2;
            } else if value > second_best {
                second_best = value;
                second_best_index = i-2;
            }
        }
        clusters.push((cell.to_string(), losses));
        //println!("cell {} index {} clusters {} {}",cell, cell_clusters.len(),best_index, second_best_index);
        cell_clusters.push((best_index, second_best_index));
    }
    (cell_clusters, clusters, num_clusters+1)
}


/// Computes log(exp(a) + exp(b)) in a way that tries to avoid over/underflows.
/// log(exp(a - b) + 1) + b = log((exp(a) + exp(b)) - exp(b)) + b = log(exp(a) + exp(b))
pub fn logaddexp(a: f64, b: f64, base: f64) -> f64 {
    let c: f64;
    let d: f64;

    // powf(a-b) doesn't do the right thing if a > b.
    if a > b {
        c = b;
        d = a;
    } else {
        c = a;
        d = b;
    }

    (base.powf(c - d) + 1.0).log(base) + d
}

#[derive(Clone)]
struct Params {
    alts: String,
    refs: String,
    clusters: String,
    p_soup: f64,
    doublet_prior: f64,
    debug: FnvHashSet<usize>,
    doublet_threshold: f64,
    singlet_threshold: f64,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let alts = params.value_of("alts").unwrap();
    let refs = params.value_of("refs").unwrap();
    let clusters = params.value_of("clusters").unwrap();
    let soup = params.value_of("p_soup").unwrap_or("0.0");
    let soup = soup.to_string().parse::<f64>().unwrap();
    let doublet_prior = params.value_of("doublet_prior").unwrap_or("0.5");
    let doublet_prior = doublet_prior.to_string().parse::<f64>().unwrap();
    let debug_cells: Vec<&str> = match params.values_of("debug") {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut debug: FnvHashSet<usize> = FnvHashSet::default();
    for cell in debug_cells {
        let cell = cell.to_string().parse::<usize>().unwrap();
        debug.insert(cell);
    }
    let doublet_threshold = params.value_of("doublet_threshold").unwrap_or("0.9");
    let doublet_threshold = doublet_threshold.to_string().parse::<f64>().unwrap();
    let singlet_threshold = params.value_of("singlet_threshold").unwrap_or("0.9");
    let singlet_threshold = singlet_threshold.to_string().parse::<f64>().unwrap();

    Params{
        alts: alts.to_string(),
        refs: refs.to_string(),
        clusters: clusters.to_string(),
        p_soup: soup,
        doublet_prior: doublet_prior,
        debug: debug,
        doublet_threshold: doublet_threshold,
        singlet_threshold: singlet_threshold,
    }
}
