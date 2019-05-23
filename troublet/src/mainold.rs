extern crate fnv;
use fnv::FnvHashMap;
use std::env;
use std::io::BufReader;
use std::io::BufRead;
use std::fs::File;

static DEBUG: bool = false;
static DEBUG_CELL: &str = "ACGGGTCGTTCCACGG-1";

fn main() {
    let args: Vec<String> = env::args().collect();
    let calls_file = &args[1];
    let cluster_file = &args[2];
    let cluster_genotype_file = &args[3];
    let p_soup = args[4].to_string().parse::<f64>().unwrap();
    let data: AlleleData = load_data(calls_file.to_string(), cluster_file.to_string(), cluster_genotype_file.to_string(), p_soup);
    let p_singlet: f64 = 0.98;
    assign_cells(data, p_singlet); 
}

fn assign_cells(data: AlleleData, p_singlet: f64) {
    let log_priors: FnvHashMap<(i32, i32),f64> = calc_log_priors(&data.cluster_cells, p_singlet);
    let mut map: FnvHashMap<String, i32> = FnvHashMap::default();
    map.insert(DEBUG_CELL.to_string(),0);

    let cells_to_assign = match DEBUG {
        true => &map,
        false => &data.cell_cluster,
    };
    for (cell, cluster_assignment) in cells_to_assign {
        let mut log_p_datas: FnvHashMap<(i32,i32),f64> = calc_log_p_datas(cell.to_string(), &data);
        let mut log_denom: f64 = 0.0;
        for (cluster_pair, log_prob) in &log_p_datas {
            if log_denom == 0.0 {
                log_denom = log_prob + log_priors.get(cluster_pair).unwrap();
            } else {
                log_denom = logaddexp(log_denom, log_prob + log_priors.get(cluster_pair).unwrap(),(1.0f64).exp());
            }
        }
        //println!("log denom {}",log_denom);
        let mut best_c1 = -1;
        let mut best_c2 = -1;
        let mut best_posterior = -1.0;
        let mut posterior_sum: f64 = 0.0;
        for (cluster_pair, log_prob) in &log_p_datas {
            //println!("{},{},{},{},{},{}", cluster_pair.0, cluster_pair.1, log_prob, log_priors.get(cluster_pair).unwrap(), (log_prob + log_priors.get(cluster_pair).unwrap() - log_denom), (log_prob + log_priors.get(cluster_pair).unwrap() - log_denom).exp());
            let posterior = (log_prob + log_priors.get(cluster_pair).unwrap() - log_denom).exp();
            posterior_sum += posterior;
            if posterior > best_posterior {
                best_c1 = cluster_pair.0;
                best_c2 = cluster_pair.1;
                best_posterior = posterior;
            }
        }
        assert!((posterior_sum - 1.0).abs() < 0.001, "posterior sum must approximately equal 1.0 but is {}", posterior_sum); // make sure our posteriors sum to 1.0
        if !DEBUG {
            println!("{}\t{}\t{}\t{}\t{}\t{}", cell, best_c1, best_c2, best_c1 != best_c2, best_posterior, cluster_assignment);
        }
    }
}

fn calc_log_p_datas(cell: String, data: &AlleleData) -> FnvHashMap<(i32, i32), f64> {
    let mut log_p_datas: FnvHashMap<(i32, i32), f64> = FnvHashMap::default();
    for (cluster1, _) in &data.cluster_cells {
        for (cluster2, _) in &data.cluster_cells {
            if cluster1 > cluster2 {
                continue
            }
            log_p_datas.insert((*cluster1, *cluster2), calc_log_p_data(&cell, &cluster1, &cluster2, &data));   
        }
    }
    log_p_datas
}

fn calc_log_p_data(cell: &String, c1: &i32, c2: &i32, data: &AlleleData) -> f64{
    let mut log_prob = 0.0;
    for allele in &data.alleles_list {
        //let c1_count = match data.cluster_allele_counts.get(c1).unwrap().get(&allele) {
        //    Some(x) => *x as f64,
        //    None => 0.0,
        //};
        //let c2_count = match data.cluster_allele_counts.get(c2).unwrap().get(&allele) {
        //    Some(x) => *x as f64,
        //    None => 0.0,
        //};
        //let comb_cluster_rate = (c1_count + c2_count + psuedocount)/
        //    ((data.cluster_cells.get(c1).unwrap() + data.cluster_cells.get(c2).unwrap()) as f64);
        //log_prob += match data.cells_to_alleles.get(cell).unwrap().get(&allele) {
        //    Some(_) => comb_cluster_rate.ln(),
        //    None => (1.0-comb_cluster_rate).ln(),
        //};
        if !data.cells_to_alleles.get(cell).unwrap().contains_key(&allele) &&
             !data.cells_to_alleles.get(cell).unwrap().contains_key(&Allele{chrom: allele.chrom, pos: allele.pos, is_ref: !allele.is_ref}){
            continue
        } 
        let log_rate = data.cluster_pair_allele_log_frequencies.get(&(*c1,*c2)).unwrap().get(allele).unwrap();
        log_prob += match data.cells_to_alleles.get(cell).unwrap().get(&allele) {
            Some(_) => log_rate.0,
            None => log_rate.1,
        };
        //if DEBUG {
        //    println!("calc_log_prob allele {} {} {}, clusters {} {} with rates c1 = {}, c2 = {}, (c1,c2) = {} and cell {} ", 
        //        allele.chrom, allele.pos, match allele.is_ref {true => "ref", false => "alt"},
        //        c1,c2, 
        //        data.cluster_pair_allele_log_frequencies.get(&(*c1,*c1)).unwrap().get(allele).unwrap().0.exp(),
        //        data.cluster_pair_allele_log_frequencies.get(&(*c2,*c2)).unwrap().get(allele).unwrap().0.exp(),
        //        data.cluster_pair_allele_log_frequencies.get(&(*c1,*c2)).unwrap().get(allele).unwrap().0.exp(),
        //        match data.cells_to_alleles.get(cell).unwrap().get(&allele) { Some(_) => "allele present", None => "allele absent"}
        //        );
        //}
        
    }
    log_prob
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

fn calc_log_priors(cluster_cells: &FnvHashMap<i32, i32>, p_singlet: f64) -> FnvHashMap<(i32,i32),f64> {
    let mut priors: FnvHashMap<(i32, i32),f64> = FnvHashMap::default();
    let mut total_cells = 0;
    for (_, cells) in cluster_cells {
        total_cells += cells;
    }
    for (cluster1, cells1) in cluster_cells {
        for (cluster2, cells2) in cluster_cells {
            if cluster1 > cluster2 {
                continue
            }
            priors.insert((*cluster1, *cluster2), calc_log_prior(*cluster1, *cluster2, *cells1, *cells2, p_singlet, total_cells));
            //priors.insert((*cluster1, *cluster2), calc_log_prior_uninformed(cluster_cells.len() as i32));
        }
    }
    priors
}

fn calc_log_prior_uninformed(clusters: i32) -> f64{
    let prior = 1.0/((clusters*(clusters-1)/2 + clusters) as f64);
    prior.ln()
}

fn calc_log_prior_uninformed_singlets(clusters: i32) ->f64 {
    (1.0/(clusters as f64)).ln()
}

fn calc_log_prior(c1: i32, c2: i32, cells1: i32, cells2: i32, p_singlet: f64, total_cells: i32) -> f64{
    if c1 == c2 {
        let singlet = p_singlet * (cells1 as f64) / (total_cells as f64);
        let doublet = (1.0-p_singlet) * (cells1 as f64) / (total_cells as f64) * (cells2 as f64) / (total_cells as f64);
        return (singlet + doublet).ln();
    }
    let doublet = (1.0-p_singlet) * (cells1 as f64) / (total_cells as f64) * (cells2 as f64) / (total_cells as f64);
    doublet.ln()
}

fn load_data(cell_variants_filename: String, cluster_filename: String, cluster_genotype_filename: String, p_soup: f64) -> AlleleData {
    let cell_variants_file = File::open(cell_variants_filename).expect("file not found");    
    let cell_variants_file_reader = BufReader::new(&cell_variants_file);
    let mut chromosomes: FnvHashMap<String, i32> = FnvHashMap::default();
    let mut chromosome_index: i32 = 0;

    let mut theta: FnvHashMap<Allele,Vec<f64>> = FnvHashMap::default();
    let cluster_genotype_file = File::open(cluster_genotype_filename).expect("file not found");
    let cluster_genotype_reader = BufReader::new(&cluster_genotype_file);
    for line in cluster_genotype_reader.lines() {
        let temp = line.unwrap();
        let tokens: Vec<&str> = temp.trim().split("\t").collect();
        let chrom = tokens[0];
        if let None = chromosomes.get(chrom) {
            chromosome_index += 1;
            chromosomes.insert(chrom.to_string(),chromosome_index);
        }
        let chrom = chromosomes.get(tokens[0]).unwrap();
        let pos = tokens[1].to_string().parse::<i32>().unwrap();
        let ref_allele = Allele{chrom:*chrom, pos:pos, is_ref:true};
        let alt_allele = Allele{chrom:*chrom, pos:pos, is_ref:false};
        theta.insert(ref_allele, vec![]);
        theta.insert(alt_allele, vec![]);
        for token in tokens[4..].iter() {
            let sub_tokens: Vec<&str> = token.split(",").collect();
            for (i, sub_token) in sub_tokens.iter().enumerate() {
                let mut probability = sub_token.to_string().parse::<f64>().unwrap();
                probability = probability.max(0.05).min(0.95); // lets not be overly confident in either direction
                
                if i == 0 {
                    theta.get_mut(&ref_allele).unwrap().push(probability);
                } else {
                    theta.get_mut(&alt_allele).unwrap().push(probability);
                }
            }   
        }
    } 
    let mut cells_to_alleles: FnvHashMap<String, FnvHashMap<Allele, bool>> = FnvHashMap::default();
    let mut alleles_list: Vec<Allele> = Vec::new();

    for line in cell_variants_file_reader.lines() {
        let linestring = line.unwrap();
        let tokens: Vec<&str> = linestring.trim().split("\t").collect();
        let chrom = tokens[0];
        let pos = tokens[1].to_string().parse::<i32>().unwrap();
        let cells_list: Vec<&str> = tokens[7].split(";").collect();
        let ref_cells: Vec<&str> = cells_list[0].split(",").collect();
        let alt_cells: Vec<&str> = cells_list[1].split(",").collect();
        let alt_allele = Allele {chrom: *chromosomes.get(chrom).unwrap(), pos: pos, is_ref: false};
        let ref_allele = Allele {chrom: *chromosomes.get(chrom).unwrap(), pos: pos, is_ref: true};
        alleles_list.push(ref_allele);
        alleles_list.push(alt_allele);
        for alt_cell in &alt_cells {
            
            if DEBUG && alt_cell.to_string() == DEBUG_CELL {
                println!("{} has {} {} {} theta {} {} {} {} {} {}",DEBUG_CELL, alt_allele.chrom, alt_allele.pos, "alt",
                    theta.get(&alt_allele).unwrap().get(0).unwrap(),
                    theta.get(&alt_allele).unwrap().get(1).unwrap(),
                    theta.get(&alt_allele).unwrap().get(2).unwrap(),
                    theta.get(&alt_allele).unwrap().get(3).unwrap(),
                    theta.get(&alt_allele).unwrap().get(4).unwrap(),
                    theta.get(&alt_allele).unwrap().get(5).unwrap()
                    );
            }
            let alleles = cells_to_alleles.entry(alt_cell.to_string()).or_insert(FnvHashMap::default());
            alleles.insert(alt_allele, true);
        }
        for ref_cell in &ref_cells {
            if DEBUG && ref_cell.to_string() == DEBUG_CELL {
                println!("{} has {} {} {} theta {} {} {} {} {} {}",DEBUG_CELL, ref_allele.chrom, ref_allele.pos, "ref",
                    theta.get(&ref_allele).unwrap().get(0).unwrap(),
                    theta.get(&ref_allele).unwrap().get(1).unwrap(),
                    theta.get(&ref_allele).unwrap().get(2).unwrap(),
                    theta.get(&ref_allele).unwrap().get(3).unwrap(),
                    theta.get(&ref_allele).unwrap().get(4).unwrap(),
                    theta.get(&ref_allele).unwrap().get(5).unwrap()
                
                );
            }
            let alleles = cells_to_alleles.entry(ref_cell.to_string()).or_insert(FnvHashMap::default());
            alleles.insert(ref_allele, true);
        }
    }

    let mut cell_cluster: FnvHashMap<String, i32> = FnvHashMap::default();    
    
    let cluster_file = File::open(cluster_filename).expect("file not found");    
    let cluster_file_reader = BufReader::new(&cluster_file);
    for line in cluster_file_reader.lines() {
        let temp = line.unwrap();
        let tokens: Vec<&str> = temp.trim().split(",").collect();
        let cell = tokens[0];
        let cluster = tokens[1].to_string().parse::<i32>().unwrap();
        cell_cluster.insert(cell.to_string(), cluster);
    }

   
    let mut loci_expression: FnvHashMap<(i32, i32), f64> = FnvHashMap::default();
    //let total_cells:f64 = cell_cluster.len() as f64;
    let mut cluster_allele_counts: FnvHashMap<i32, FnvHashMap<Allele, i32>> = FnvHashMap::default();
    let mut cluster_cells: FnvHashMap<i32,i32> = FnvHashMap::default();
    let mut total_allele_counts: FnvHashMap<Allele,i32> = FnvHashMap::default();
    let total_cells = cell_cluster.len();
    for (cell, alleles) in &cells_to_alleles {
        let cluster = *cell_cluster.get(cell).unwrap();
        let count = cluster_cells.entry(cluster).or_insert(0);
        *count += 1;
        cluster_allele_counts.entry(cluster).or_insert(FnvHashMap::default());
        for (allele, _) in alleles {
            let count = cluster_allele_counts.get_mut(&cluster).unwrap().entry(*allele).or_insert(0);
            *count += 1;
            let count = total_allele_counts.entry(*allele).or_insert(0);
            *count += 1;
            let count = loci_expression.entry((allele.chrom, allele.pos)).or_insert(0.0f64);
            *count += 1.0f64/(total_cells as f64);
        }
    }
    let mut cluster_pair_allele_log_frequencies: FnvHashMap<(i32, i32), FnvHashMap<Allele, (f64, f64)>> = FnvHashMap::default();
    for (c1, c1_allele_counts) in &cluster_allele_counts {
        for (c2, c2_allele_counts) in &cluster_allele_counts {
            cluster_pair_allele_log_frequencies.insert((*c1,*c2), FnvHashMap::default());
            for allele in &alleles_list {
                let c1_counts = match c1_allele_counts.get(allele) {
                    Some(x) => *x as f64,
                    None => 0.0f64, 
                };
                let c2_counts = match c2_allele_counts.get(allele) {
                    Some(x) => *x as f64,
                    None => 0.0f64,
                };
                let theta_prob1 = theta.get(allele).unwrap()[*c1 as usize];
                let theta_prob2 = theta.get(allele).unwrap()[*c2 as usize];
                
                let p1 = c1_counts/(*cluster_cells.get(c1).unwrap() as f64) * theta_prob1;
                let p2 = c2_counts/(*cluster_cells.get(c2).unwrap() as f64) * theta_prob2;
                //let p1 = loci_expression.get(&(allele.chrom, allele.pos)).unwrap() * theta_prob1;
                //let p2 = loci_expression.get(&(allele.chrom, allele.pos)).unwrap() * theta_prob2;
                //let p1 = 0.5*theta_prob1;
                //let p2 = 0.5*theta_prob2;

                let allele_rate = (*total_allele_counts.get(allele).unwrap() as f64)/(total_cells as f64);
                let rate = match c1 == c2 {
                    true => p_soup*allele_rate + (1.0-p_soup)*(1.0-(1.0-theta_prob1)*(1.0-theta_prob2)),//p_soup*allele_rate + (1.0 - p_soup) * p1,
                    //false => p_soup*((1.0 - theta_prob1) * (1.0 - theta_prob2)) + (1.0 - p_soup) * (1.0 - (1.0 - theta_prob1) * (1.0 - theta_prob2)),
                    false => p_soup*allele_rate + (1.0-p_soup)*(1.0-(1.0-theta_prob1)*(1.0-theta_prob2)),
                };
                //if DEBUG {
                //    if cells_to_alleles.get(DEBUG_CELL).unwrap().contains_key(allele) {
                //        println!("cluster rates c1 = {}, c2 = {}, allele = chrom {}, pos {}, is_ref {}, has counts {} {} and sizes {} {} for rate {}",
                //            c1, c2, allele.chrom, allele.pos, match allele.is_ref {true => "ref", false => "alt",},
                //            c1_counts, c2_counts,
                //            *cluster_cells.get(c1).unwrap(), *cluster_cells.get(c2).unwrap(),
                //            rate,
                //        );
                //    }
                //}             
                //(c1_counts + c2_counts + psuedocount)/(*cluster_cells.get(c1).unwrap() as f64 + *cluster_cells.get(c2).unwrap() as f64);
                cluster_pair_allele_log_frequencies.get_mut(&(*c1,*c2)).unwrap().insert(*allele, (rate.ln(),(1.0-rate).ln()));
            }
        }
    }

    AlleleData{
        cells_to_alleles: cells_to_alleles,
        alleles_list: alleles_list,
        cell_cluster: cell_cluster,
        //cluster_allele_counts: cluster_allele_counts,
        cluster_cells: cluster_cells,   
        cluster_pair_allele_log_frequencies: cluster_pair_allele_log_frequencies,
    }
}

struct AlleleData {
    cells_to_alleles: FnvHashMap<String, FnvHashMap<Allele, bool>>,
    alleles_list: Vec<Allele>,
    cell_cluster: FnvHashMap<String, i32>,
    //cluster_allele_counts: HashMap<i32, HashMap<Allele,i32>>,
    cluster_cells: FnvHashMap<i32, i32>,
    cluster_pair_allele_log_frequencies: FnvHashMap<(i32, i32), FnvHashMap<Allele, (f64, f64)>>,
}

#[derive(PartialEq, Eq, Hash, Clone, Copy)]
struct Allele {
    chrom: i32,
    pos: i32,
    is_ref: bool,
}
