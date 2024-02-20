#[macro_use]
extern crate clap;
extern crate hashbrown;
extern crate rand;
extern crate statrs;
extern crate itertools;
extern crate rayon;
extern crate vcf;
extern crate flate2;



//Reza
extern crate csv;
use csv::Writer;
use std::error::Error;



use std::fs::OpenOptions;
use std::io::prelude::*;


use flate2::read::GzDecoder;
use flate2::read::MultiGzDecoder;
use vcf::*;


use rayon::prelude::*;

use rand::Rng;
use rand::rngs::StdRng;
use rand::SeedableRng;

use clap::App;
use std::f32;

use std::ffi::OsStr;
use std::io::Read;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use hashbrown::{HashMap,HashSet};
use itertools::izip;




fn main() {
    let params = load_params();
    let mut cell_barcodes = load_barcodes(&params); 
    io::stdout().flush().unwrap();
    let (loci_used, total_cells, cell_data, index_to_locus, locus_to_index, list_of_loci_used, locus_cell_counts, ref_alt_counts_per_locus) = load_cell_data(&params);

    let (mut cell_hashing_cluster_centers,mut cell_id_to_label_map, 
            barcode_to_cell_id_and_hash, cell_id_to_barcode_and_hash) = initialize_cluster_centers_cell_hashing(
        &list_of_loci_used, 
        &cell_data, 
        cell_barcodes.clone(), 
        &params
    );

    cellector(ref_alt_counts_per_locus,index_to_locus, &list_of_loci_used, locus_cell_counts, loci_used, cell_data, &params, cell_barcodes, locus_to_index,cell_id_to_barcode_and_hash);

    // souporcell_main(&list_of_loci_used,locus_cell_counts,cell_hashing_cluster_centers,loci_used, cell_data, &params, cell_barcodes, locus_to_index, cell_id_to_barcode_and_hash);
}
   






struct ThreadData {
    best_log_probabilities: Vec<Vec<f32>>,
    best_total_log_probability: f32,
    rng: StdRng,
    solves_per_thread: usize,
    thread_num: usize,
}

impl ThreadData {
    fn from_seed(seed: [u8; 32], solves_per_thread: usize, thread_num: usize) -> ThreadData {
        ThreadData {
            best_log_probabilities: Vec::new(),
            best_total_log_probability: f32::NEG_INFINITY,
            rng: SeedableRng::from_seed(seed),
            solves_per_thread: solves_per_thread,
            thread_num: thread_num,
        }
    }
}




fn cellector( ref_alt_counts_per_locus:HashMap<usize, [u32; 2]>,index_to_locus: Vec<usize> ,list_of_loci_used: &HashSet<usize> ,loci_cell_count :HashMap<usize, usize> , loci: usize, cell_data: Vec<CellData>, params: &Params, barcodes: Vec<String>, locus_to_index: HashMap<usize, usize>, cell_id_to_barcode_and_hash: Vec<(String, String)>) {

    let seed = [params.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);

    let mut beta_cluster_centers = init_beta_binomial_cluster_centers(ref_alt_counts_per_locus,loci, &cell_data, params, &mut rng, index_to_locus.clone());

    let (log_loss, log_probabilities) = beta_EM(index_to_locus,&mut beta_cluster_centers.clone() , &cell_data, params, loci);
    // println!("log_loss = {}", log_loss);    
    // println!("log_probabilities = {:?}", log_probabilities);

    let mut Soup_Assignments = OpenOptions::new()
    .write(true)
    .truncate(true)
    .create(true)
    .open("assignments.tsv")
    .unwrap();

    // write header to cells_file
    writeln!(Soup_Assignments, "Barcode\tCell\tSoup").unwrap();

    let mut cell_id = 0;

    let threshold = 0.55;

    for (bc, log_probs) in barcodes.iter().zip(log_probabilities.iter()) {
        
        let mut best = 0;
        let mut best_lp = f64::NEG_INFINITY;
        for index in 0..log_probs.len() {
            if log_probs[index] > best_lp {
                best = index;
                best_lp = log_probs[index];
            }
        }
        // let posterior = (best_lp - log_sum_exp(log_probs)).exp();
        let posterior = (best_lp - log_sum_exp_64(log_probs.clone())).exp();

        writeln!(Soup_Assignments, "{}\t{}\t{}", bc, cell_id_to_barcode_and_hash[cell_id].1, best).unwrap();

        // if posterior > threshold {
        //     print!("{}\t{}\t{}",bc, best, cell_id_to_barcode_and_hash[cell_id].1);
        //     writeln!(Soup_Assignments, "{}\t{}\t{}", bc, cell_id_to_barcode_and_hash[cell_id].1, best).unwrap();
        // }
        // else {
        //     print!("{}\t{}\tU",bc,cell_id_to_barcode_and_hash[cell_id].1);
        //     writeln!(Soup_Assignments, "{}\t{}\tU", bc, cell_id_to_barcode_and_hash[cell_id].1).unwrap();
        // }
        
        
        // print!("\n");

        cell_id += 1;


    }

    




        // if posterior > threshold {

        //     print!("{}\t{}\t{}",bc, best, cell_id_to_barcode_and_hash[cell_id].1);
        //     writeln!(Soup_Assignments, "{}\t{}\t{}", bc, cell_id_to_barcode_and_hash[cell_id].1, best).unwrap();

        // } else {

        //     print!("{}\tU\t{}",bc, cell_id_to_barcode_and_hash[cell_id].1);
        //     writeln!(Soup_Assignments, "{}\t{}\tU", bc, cell_id_to_barcode_and_hash[cell_id].1 ).unwrap();
        // }

}

fn souporcell_main(list_of_loci_used: &HashSet<usize> ,loci_cell_count :HashMap<usize, usize> ,cell_hashing_assignments:  Vec<Vec<f32>>, loci_used: usize, cell_data: Vec<CellData>, params: &Params, barcodes: Vec<String>, locus_to_index: HashMap<usize, usize>, cell_id_to_barcode_and_hash: Vec<(String, String)>) {
    let seed = [params.seed; 32];
    let mut rng: StdRng = SeedableRng::from_seed(seed);
    let mut threads: Vec<ThreadData> = Vec::new();
    let solves_per_thread = ((params.restarts as f32)/(params.threads as f32)).ceil() as usize;
    for i in 0..params.threads {
        threads.push(ThreadData::from_seed(new_seed(&mut rng), solves_per_thread, i));
    }
    
    threads.par_iter_mut().for_each(|thread_data| {
        for iteration in 0..thread_data.solves_per_thread {

            let cluster_centers: Vec<Vec<f32>> = init_cluster_centers(loci_used, &cell_data, params, &mut thread_data.rng, &locus_to_index);
            // let cluster_centers = cell_hashing_assignments.clone();
            
            let (log_loss, log_probabilities) = EM(locus_to_index.clone() ,list_of_loci_used ,loci_cell_count.clone() ,loci_used, cluster_centers, &cell_data ,params, iteration, thread_data.thread_num, cell_id_to_barcode_and_hash.clone());            



            if log_loss > thread_data.best_total_log_probability {
                thread_data.best_total_log_probability = log_loss;
                thread_data.best_log_probabilities = log_probabilities;
            }
        
         
        }
    });
    let mut best_log_probability = f32::NEG_INFINITY;
    let mut best_log_probabilities: Vec<Vec<f32>> = Vec::new();
    for thread_data in threads {
        if thread_data.best_total_log_probability > best_log_probability {
            best_log_probability = thread_data.best_total_log_probability;
            best_log_probabilities = thread_data.best_log_probabilities;
        }
    }
    eprintln!("best total log probability = {}", best_log_probability);



    let mut cell_id = 0;
    let mut num_unassigned = 0;
    let threshold = 0.99;


    let mut Soup_Assignments = OpenOptions::new()
    .write(true)
    .truncate(true)
    .create(true)
    .open("assignments.tsv")
    .unwrap();

    // write header to cells_file
    writeln!(Soup_Assignments, "Barcode\tCell\tSoup").unwrap();



    for (bc, log_probs) in barcodes.iter().zip(best_log_probabilities.iter()) {
        let mut best = 0;
        let mut best_lp = f32::NEG_INFINITY;
        for index in 0..log_probs.len() {
            if log_probs[index] > best_lp {
                best = index;
                best_lp = log_probs[index];
                
            }
        }
        let posterior = (best_lp - log_sum_exp(log_probs)).exp();




        if posterior > threshold {

            print!("{}\t{}\t{}",bc, best, cell_id_to_barcode_and_hash[cell_id].1);
            writeln!(Soup_Assignments, "{}\t{}\t{}", bc, cell_id_to_barcode_and_hash[cell_id].1, best).unwrap();

        } else {

            print!("{}\tU\t{}",bc, cell_id_to_barcode_and_hash[cell_id].1);
            writeln!(Soup_Assignments, "{}\t{}\tU", bc, cell_id_to_barcode_and_hash[cell_id].1 ).unwrap();
        }
        for index in 0..log_probs.len() {
            print!("{}",log_probs[index]);
            if index < log_probs.len() - 1 { print!("\t"); } 
    
        } print!("\n");
    

        cell_id += 1;

    }

}



// REZA EM function for beta binomial
fn beta_EM(index_to_locus: Vec<usize>, beta_cluster_centers: &mut Vec<Vec<(f32, f32)>>, cell_data: &Vec<CellData>,params: &Params, loci: usize) -> (f64, Vec<Vec<f64>>) {

    let mut beta_loss_values: Vec<f64> = Vec::new();
    let mut total_log_loss = f32::NEG_INFINITY;
    let mut beta_loss_values: Vec<f64> = Vec::new();
    let mut iterations = 0;
    let mut last_log_loss = f64::NEG_INFINITY;
    let mut log_loss_change = 10000.0;
    let log_loss_change_limit = 0.0001*(cell_data.len() as f64);
    let mut log_beta_loss_total = 0.0;
    let mut cell_cluster_probabilities: Vec<Vec<f64>> = vec![vec![0.0; 2]; cell_data.len()];

    let mut final_log_probabilities: Vec<Vec<f64>> = Vec::new();
    for _cell in 0..cell_data.len() {
        final_log_probabilities.push(Vec::new());
    }



    while iterations < 1 {
        
        log_beta_loss_total = 0.0;
        for (celldex, cell) in cell_data.iter().enumerate() {

            beta_loss_values = log_beta_loss(cell, &beta_cluster_centers, celldex, params);
            println!("beta_loss_values = {:?}", beta_loss_values);

            log_beta_loss_total += log_sum_exp_64(beta_loss_values.clone());

            let normalize_in_log = normalize_log_probabilities2(beta_loss_values);
            final_log_probabilities[celldex] = normalize_in_log.clone();
            // println!("final logs= {:?}", final_log_probabilities[celldex]);   
        }
        

        update_alpha_and_beta(index_to_locus.clone(), beta_cluster_centers, cell_data, &mut cell_cluster_probabilities, loci, &final_log_probabilities, params);
        println!("loss = {}", log_beta_loss_total);
        log_loss_change = log_beta_loss_total - last_log_loss;
        last_log_loss = log_beta_loss_total;

        


        iterations += 1;
    }
    




    (log_beta_loss_total, final_log_probabilities)

}

//REZA normalization for beta 
fn normalize_log_probabilities(log_probs: Vec<f64>) -> Vec<f64> {
    let max_log_prob = log_probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let sum_exp_adjusted: f64 = log_probs
        .iter()
        .map(|&log_prob| (log_prob - max_log_prob).exp())
        .sum();
    let log_sum_exp_adjusted = sum_exp_adjusted.ln();

    log_probs
        .iter()
        .map(|&log_prob| log_prob - (log_sum_exp_adjusted + max_log_prob))
        .collect()
}

fn normalize_log_probabilities2(log_probs: Vec<f64>) -> Vec<f64> {
    let mut normalized_probabilities: Vec<f64> = Vec::new();
    let sum = log_sum_exp_64(log_probs.clone());
    for i in 0..log_probs.len() {
        normalized_probabilities.push((log_probs[i]-sum).exp());
    }
    normalized_probabilities
}


//REZA beta update
fn update_alpha_and_beta(
    index_to_locus: Vec<usize>,
    cluster_centers: &mut Vec<Vec<(f32, f32)>>,
    cell_data: &Vec<CellData>,
    log_probabilities: &mut Vec<Vec<f64>>, 
    loci: usize,
    cell_cluster_probabilities: &Vec<Vec<f64>>,
    params: &Params
)  {

    for cluster in 0..params.num_clusters {

        for locus_index in 0..loci {
            
            let locus = index_to_locus[locus_index as usize];
            let mut updated_alpha = 1.0 as f32;
            let mut updated_beta = 1.0 as f32;

            // println!("locus index = {}", locus_index);
            for (celldex, cell) in cell_data.iter().enumerate(){
                
                if cell.loci.contains(&locus) {
                    let mut alt_counts = 0;
                    let mut ref_counts = 0;
                    alt_counts = cell.alt_counts[cell.loci.iter().position(|&r| r == locus).unwrap()] as usize;
                    ref_counts = cell.ref_counts[cell.loci.iter().position(|&r| r == locus).unwrap()] as usize;
                    updated_alpha += cell_cluster_probabilities[celldex][cluster] as f32 * alt_counts as f32;
                    updated_beta += cell_cluster_probabilities[celldex][cluster] as f32 * ref_counts as f32;
                }

            }

            cluster_centers[cluster][locus_index] = (updated_alpha, updated_beta);

        }
    }
}




fn EM(locus_to_index: HashMap<usize, usize> ,list_of_loci_used: &HashSet<usize> ,loci_cell_count :HashMap<usize, usize>, loci: usize, mut cluster_centers: Vec<Vec<f32>>, cell_data: &Vec<CellData>, params: &Params, epoch: usize, thread_num: usize, cell_id_to_barcode_and_hash: Vec<(String, String)>) -> (f32, Vec<Vec<f32>>) {
    let mut sums: Vec<Vec<f32>> = Vec::new();
    let mut denoms: Vec<Vec<f32>> = Vec::new();


    for cluster in 0..params.num_clusters {
        sums.push(Vec::new());
        denoms.push(Vec::new());
        for index in 0..loci {
            sums[cluster].push(1.0);
            denoms[cluster].push(2.0); // HAYNES debug, no psuedocounts
        }
    }
    

    let mut file = OpenOptions::new()
    .write(true)
    .append(true)
    .create(true)
    .open("log_loss_change_with_last.txt")
    .unwrap();



    let mut log_binoms: Vec<f32> = Vec::new();
    let mut beta_loss_values: Vec<f64> = Vec::new();


    // let mut cells_file = OpenOptions::new()
    // .write(true)
    // .append(true)
    // .create(true)
    // .open("cells_contribution.tsv")
    // .unwrap();

    // // write header to cells_file
    // writeln!(cells_file, "cell_id\tlocus_index\tref_counts\talt_counts\tlog_loss\tcluster").unwrap();




    let mut log_priors: Vec<f32> = Vec::new();
    for i in 0..params.num_clusters { log_priors.push((1.0/(params.num_clusters as f32)).ln()); }


    let mut iterations = 0;

    let mut total_log_loss = f32::NEG_INFINITY;

    let mut final_log_probabilities = Vec::new();
    for _cell in 0..cell_data.len() {
        final_log_probabilities.push(Vec::new());
    }
    let log_loss_change_limit = 0.0001*(cell_data.len() as f32);
    let temp_steps = 1; // debug no deterministic annealing
    let mut last_log_loss = f32::NEG_INFINITY;


    let mut last_log_loss_binom = f32::NEG_INFINITY;
    let mut count_wrong_directions = 0;
    



    for temp_step in 0..temp_steps {
        // println!("temp_step={}", temp_step);
      
        let mut log_loss_change = 10000.0;
    
        let mut temp = 1.0;


        while log_loss_change > log_loss_change_limit && iterations < 10000 {


           
        
            let mut log_binom_loss = 0.0;
            let mut log_beta_loss_total = 0.0;
           
            let mut total_posteriors: Vec<f32> = Vec::new();
            for i in 0..params.num_clusters { total_posteriors.push(0.0); }


            // let mut cells_file = OpenOptions::new()
            // .write(true)
            // .append(true)
            // .create(true)
            // .open("cells_contribution_".to_owned() + &iterations.to_string() + ".tsv")
            // .unwrap();
            // writeln!(cells_file, "cell_id\tlocus_index\tref_counts\talt_counts\tlog_loss\tcluster").unwrap();
            reset_sums_denoms(loci, &mut sums, &mut denoms, &cluster_centers, params.num_clusters);



            for (celldex, cell) in cell_data.iter().enumerate() {

                
                // let log_binoms = binomial_loss(&mut cells_file, cell, &cluster_centers, &log_priors, celldex);
                // let log_binoms = binomial_loss(cell, &cluster_centers, &log_priors, celldex);
                log_binoms = binomial_loss(cell, &cluster_centers, &log_priors, celldex);
                // let log_binoms = sum_of_squares_loss(cell, &cluster_centers, log_priors[0], celldex);

               


                log_binom_loss += log_sum_exp(&log_binoms);
                temp = (cell.total_alleles/(20.0 * 2.0f32.powf(temp_step as f32))).max(1.0);
                if temp_step == temp_steps - 1 { temp = 1.0; }

                //if temp_step > 0 { temp = 1.0; }

                let mut probabilities = normalize_in_log_with_temp(&log_binoms, temp);
                
                for (i, prob) in probabilities.iter().enumerate() {
                    total_posteriors[i] += prob; // this is for later updating priors
                }

                update_centers_average(&mut sums, &mut denoms, cell, &probabilities);
                final_log_probabilities[celldex] = log_binoms.clone();//log_probabilities;                


            }
            

            total_log_loss = log_binom_loss;
            log_loss_change = log_binom_loss - last_log_loss;//log_loss - last_log_loss;
            last_log_loss = log_binom_loss;//log_loss;


            update_final(loci, &sums, &denoms, &mut cluster_centers);
            let mut sum_posteriors: f32 = 0.0;
            for postsum in &total_posteriors { sum_posteriors += postsum; }
            for i in 0..params.num_clusters {
                log_priors[i] = (total_posteriors[i]/sum_posteriors).ln();
            }

            // wtr.flush();
            iterations += 1;

            
        }
        println!(" log_loss = {}", total_log_loss);

        

        
    }

    (total_log_loss, final_log_probabilities)
}




//REZA beta_loss function
fn log_beta_loss(cell_data: &CellData, cluster_centers: &Vec<Vec<(f32, f32)>>, _cellnum: usize, params: &Params) -> Vec<f64> {

    let mut probabilities: Vec<f64> = Vec::new();
    for (cluster, centers) in cluster_centers.iter().enumerate() {
        let mut prior: f64 = (1.0/params.num_clusters as f64).ln();
        probabilities.push(prior);
        for (locus_index, locus) in cell_data.loci.iter().enumerate() {

            let alt_counts = cell_data.alt_counts[locus_index] as f32;
            let ref_counts = cell_data.ref_counts[locus_index] as f32;
            let center = &centers[locus_index];
            println!("center = {:?} for cluster : {}", center, cluster);
            let alpha = center.0;
            let beta = center.1;
            let coefficient = cell_data.log_binomial_coefficient[locus_index] as f64;
            let log_prob_num = log_beta_calc(alt_counts + alpha, ref_counts + beta);
            let log_prob_denom = log_beta_calc(alpha, beta);
            let log_prob = ( (coefficient + log_prob_num) - log_prob_denom );
            probabilities[cluster] += log_prob;

        }

    }
    // println!("probabilities = {:?}", probabilities);
    probabilities
}


fn log_beta_calc(alpha: f32, beta: f32) -> f64 {

    let log_gamma_alpha = statrs::function::gamma::ln_gamma(alpha as f64);
    let log_gamma_beta = statrs::function::gamma::ln_gamma(beta as f64);
    let log_gamma_alpha_beta = statrs::function::gamma::ln_gamma((alpha + beta) as f64);

    log_gamma_alpha + log_gamma_beta - log_gamma_alpha_beta
}



// fn binomial_loss(cells_file: &mut File, cell_data: &CellData, cluster_centers: &Vec<Vec<f32>>, log_priors: &Vec<f32>, cellnum: usize) -> Vec<f32> {

fn binomial_loss(cell_data: &CellData, cluster_centers: &Vec<Vec<f32>>, log_priors: &Vec<f32>, cellnum: usize) -> Vec<f32> {

    let mut log_probabilities: Vec<f32> = Vec::new();
    let mut sum = 0.0;
    
    for (cluster, center) in cluster_centers.iter().enumerate() {
        log_probabilities.push(log_priors[cluster]);
        for (locus_index, locus) in cell_data.loci.iter().enumerate() {
            let mut prob = cell_data.log_binomial_coefficient[locus_index] + 
            (cell_data.alt_counts[locus_index] as f32) * center[*locus].ln() + 
            (cell_data.ref_counts[locus_index] as f32) * (1.0 - center[*locus]).ln();
            log_probabilities[cluster] += prob;
            // writeln!(cells_file, "{}\t{}\t{}\t{}\t{}\t{}",cellnum, cell_data.loci[locus_index], cell_data.ref_counts[locus_index], cell_data.alt_counts[locus_index], prob, cluster).unwrap();

        }
  
        sum += log_probabilities[cluster];
    }

    log_probabilities
}






fn sum_of_squares_loss(cell_data: &CellData, cluster_centers: &Vec<Vec<f32>>, log_prior: f32, cellnum: usize) -> Vec<f32> {
    let mut log_probabilities: Vec<f32> = Vec::new();
    for (cluster, center) in cluster_centers.iter().enumerate() {
        log_probabilities.push(log_prior);
        for (locus_index, locus) in cell_data.loci.iter().enumerate() {
            log_probabilities[cluster] -= (cell_data.allele_fractions[locus_index] - center[*locus]).powf(2.0);
        }
    }
    log_probabilities 
}



//REZA
fn log_sum_exp_64(p: Vec<f64>) -> f64{
    let max_p: f64 = p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let sum_rst: f64 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln() 
}


fn log_sum_exp(p: &Vec<f32>) -> f32{
    let max_p: f32 = p.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    let sum_rst: f32 = p.iter().map(|x| (x - max_p).exp()).sum();
    max_p + sum_rst.ln()
}

fn normalize_in_log(log_probs: &Vec<f32>) -> Vec<f32> { // takes in a log_probability vector and converts it to a normalized probability
    let mut normalized_probabilities: Vec<f32> = Vec::new();
    let sum = log_sum_exp(log_probs);
    for i in 0..log_probs.len() {
        normalized_probabilities.push((log_probs[i]-sum).exp());
    }
    normalized_probabilities
}


fn normalize_in_log_with_temp(log_probs: &Vec<f32>, temp: f32) -> Vec<f32> {
    let mut normalized_probabilities: Vec<f32> = Vec::new();
    let mut new_log_probs: Vec<f32> = Vec::new();
    for log_prob in log_probs {
        new_log_probs.push(log_prob/temp);
    }
    let sum = log_sum_exp(&new_log_probs);
    for i in 0..log_probs.len() {
        normalized_probabilities.push((new_log_probs[i]-sum).exp());
    }
    normalized_probabilities 
}

fn update_final(loci: usize, sums: &Vec<Vec<f32>>, denoms: &Vec<Vec<f32>>, cluster_centers: &mut Vec<Vec<f32>>) {
    for locus in 0..loci {
        for cluster in 0..sums.len() {
            let update = sums[cluster][locus]/denoms[cluster][locus];
            cluster_centers[cluster][locus] = update;//.min(0.99).max(0.01);//max(0.0001, min(0.9999, update));
        }
    }
}

fn reset_sums_denoms(loci: usize, sums: &mut Vec<Vec<f32>>, 
    denoms: &mut Vec<Vec<f32>>, cluster_centers: &Vec<Vec<f32>>, num_clusters: usize) {
    for cluster in 0..num_clusters {
        for index in 0..loci {
            sums[cluster][index] = 1.0;
            denoms[cluster][index] = 2.0;
        }
    }
}


fn update_centers_flat(sums: &mut Vec<Vec<f32>>, denoms: &mut Vec<Vec<f32>>, cell: &CellData, probabilities: &Vec<f32>) {
    for locus in 0..cell.loci.len() {
        for (cluster, probability) in probabilities.iter().enumerate() {
            sums[cluster][cell.loci[locus]] += probabilities[cluster] * cell.allele_fractions[locus];
            denoms[cluster][cell.loci[locus]] += probabilities[cluster];
        }
    }
}

fn update_centers_average(sums: &mut Vec<Vec<f32>>, denoms: &mut Vec<Vec<f32>>, cell: &CellData, probabilities: &Vec<f32>) {
    for locus in 0..cell.loci.len() {
        for (cluster, probability) in probabilities.iter().enumerate() {
            sums[cluster][cell.loci[locus]] += probabilities[cluster] * (cell.alt_counts[locus] as f32);
            denoms[cluster][cell.loci[locus]] += probabilities[cluster] * ((cell.alt_counts[locus] + cell.ref_counts[locus]) as f32);
        }
    }
}




//REZA beta binomial cluster centers initialization
fn init_beta_binomial_cluster_centers(ref_alt_counts_per_locus:HashMap<usize, [u32; 2]>, loci: usize, cell_data: &Vec<CellData>, params: &Params, rng: &mut StdRng, index_to_locus: Vec<usize>) -> Vec<Vec<(f32, f32)>> {

    let mut centers: Vec<Vec<(f32, f32)>> = Vec::new();

    for cluster in 0..params.num_clusters {
        centers.push(Vec::new());
        for _ in 0..loci {
            // Generate two random f32 values for each tuple, ensuring they are within your specified range
            let center_tuple = (
                rng.gen::<f32>().min(0.0001).max(0.0002),
                rng.gen::<f32>().min(0.0001).max(0.0002)
            );
            centers[cluster].push(center_tuple);
        }
    }



    for locus_index in 0..loci {



        //Counting ref and alt per loci

        let mut ref_counts = 0.0 as f32;
        let mut alt_counts = 0.0 as f32;
        
        let locus = index_to_locus[locus_index as usize];

        let counts = ref_alt_counts_per_locus.get(&locus).unwrap(); 
        let ref_count = counts[0] as f32; 
        let alt_count = counts[1] as f32;

        ref_counts = ref_count;
        alt_counts = alt_count;
        // println!("alt_counts = {}\t ref_counts = {}\t\t locus: {}", alt_counts, ref_counts, locus);

        let x = rng.gen_range(0.0, alt_counts);
        let y = rng.gen_range(0.0, ref_counts);

        let alpha_1 = 1.0 + x;
        let beta_1 =  1.0 + y;

        let mut alpha_2 = 1.0 + (alt_counts - x);
        let mut beta_2 = 1.0 + (ref_counts - y);

        centers[0][locus_index] = (alpha_1, beta_1);
        centers[1][locus_index] = (alpha_2, beta_2);
        // centers[0].push((alpha_1, beta_1));
        // centers[1].push((alpha_2, beta_2));

        // println!("alpha_1 = {}\t beta_1 = {}\t\t alpha_2 = {}\t beta_2 = {}", alpha_1, beta_1, alpha_2, beta_2);


    }

    
    centers

}




fn init_cluster_centers(loci_used: usize, cell_data: &Vec<CellData>, params: &Params, rng: &mut StdRng, locus_to_index: &HashMap<usize, usize>) -> Vec<Vec<f32>> {
    if let Some(known_genotypes) = &params.known_genotypes {
        return init_cluster_centers_known_genotypes(loci_used, params, rng, locus_to_index);
    } else if let Some(assigned_cells) = &params.known_cell_assignments {
        return init_cluster_centers_known_cells(loci_used, &cell_data, params, rng);
    } else {
        match params.initialization_strategy {
            ClusterInit::KmeansPP => init_cluster_centers_kmeans_pp(loci_used, &cell_data, params, rng),
            ClusterInit::RandomUniform => init_cluster_centers_uniform(loci_used, params, rng),
            ClusterInit::RandomAssignment => init_cluster_centers_random_assignment(loci_used, &cell_data, params, rng),
            ClusterInit::MiddleVariance => init_cluster_centers_middle_variance(loci_used, &cell_data, params, rng),
        }
    }
}

pub fn reader(filename: &str) -> Box<dyn BufRead> {
    let path = Path::new(filename);
    let file = match File::open(&path) {
        Err(why) => panic!("couldn't open file {}", filename),
        Ok(file) => file,
    };
    if path.extension() == Some(OsStr::new("gz")) {
        Box::new(BufReader::with_capacity(128 * 1024, MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::with_capacity(128 * 1024, file))
    }
}


fn init_cluster_centers_known_genotypes(loci: usize, params: &Params, rng: &mut StdRng, locus_to_index: &HashMap<usize, usize>) -> Vec<Vec<f32>> {
    let mut centers: Vec<Vec<f32>> = Vec::new();
    for cluster in 0..params.num_clusters {
        centers.push(Vec::new());
        for _ in 0..loci {
            centers[cluster].push(0.5);
        }
    }
    let mut vcf_reader = VCFReader::new(reader(params.known_genotypes.as_ref().unwrap())).unwrap();
    let mut locus_id: usize = 0;

    for record in vcf_reader {
        let record = record.unwrap();
        if let Some(loci_index) = locus_to_index.get(&locus_id) {
            if params.known_genotypes_sample_names.len() > 0 {
                for (sample_index, sample) in params.known_genotypes_sample_names.iter().enumerate() {
                    let gt = record.call[sample]["GT"][0].to_string();
                    // complicated way of getting the haplotype to numbers
                    let hap0 = gt.chars().nth(0).unwrap().to_string();
                    if hap0 == "." { continue; }
                    let hap0 = hap0.parse::<u32>().unwrap().min(1);
                    let hap1 = gt.chars().nth(2).unwrap().to_string().parse::<u32>().unwrap().min(1);
                    centers[sample_index][*loci_index] = (((hap0 + hap1) as f32)/2.0).min(0.99).max(0.01);
                }
            } else { assert!(false, "currently requiring known_genotypes_sample_names if known_genotypes set"); }
        }
        locus_id += 1;
    }
    centers
}

fn init_cluster_centers_known_cells(loci: usize, cell_data: &Vec<CellData>, params: &Params, rng: &mut StdRng) -> Vec<Vec<f32>> {
    assert!(false, "known cell assignments not yet implemented");
    Vec::new()
}

fn init_cluster_centers_kmeans_pp(loci: usize, cell_data: &Vec<CellData>, params: &Params, rng: &mut StdRng) -> Vec<Vec<f32>> {
    assert!(false, "kmeans++ not yet implemented");
    Vec::new()
}






fn init_cluster_centers_uniform(loci: usize, params: &Params, rng: &mut StdRng) -> Vec<Vec<f32>> {
    let mut centers: Vec<Vec<f32>> = Vec::new();
    for cluster in 0..params.num_clusters {
        centers.push(Vec::new());
        for _ in 0..loci {
            centers[cluster].push(rng.gen::<f32>().min(0.9999).max(0.0001));
        }
    }
    centers
}

fn init_cluster_centers_random_assignment(loci: usize, cell_data: &Vec<CellData>, params: &Params, rng: &mut StdRng) -> Vec<Vec<f32>> {
    let mut sums: Vec<Vec<f32>> = Vec::new();
    let mut denoms: Vec<Vec<f32>> = Vec::new();
    for cluster in 0..params.num_clusters {
        sums.push(Vec::new());
        denoms.push(Vec::new());
        for _ in 0..loci {
            sums[cluster].push(rng.gen::<f32>()*0.01);
            denoms[cluster].push(0.01);
        }
    }
    for cell in cell_data {
        let cluster = rng.gen_range(0,params.num_clusters);
        for locus in 0..cell.loci.len() {
            let alt_c = cell.alt_counts[locus] as f32;
            let total = alt_c + (cell.ref_counts[locus] as f32);
            let locus_index = cell.loci[locus];
            sums[cluster][locus_index] += alt_c;
            denoms[cluster][locus_index] += total;
        }
    }
    for cluster in 0..params.num_clusters {
        for locus in 0..loci {
            sums[cluster][locus] = sums[cluster][locus]/denoms[cluster][locus] + (rng.gen::<f32>()/2.0 - 0.25);
            sums[cluster][locus] = sums[cluster][locus].min(0.9999).max(0.0001);
        }
    }
    let centers = sums;
    
    centers
}

fn init_cluster_centers_middle_variance(loci: usize, cell_data: &Vec<CellData>, params: &Params, rng: &mut StdRng) -> Vec<Vec<f32>> {
    assert!(false, "middle variance not yet implemented");
    Vec::new()
}

fn load_cell_data(params: &Params) -> (usize, usize, Vec<CellData>, Vec<usize>, HashMap<usize, usize>, HashSet<usize>, HashMap<usize, usize>, HashMap<usize, [u32; 2]>) {
    let alt_reader = File::open(params.alt_mtx.to_string()).expect("cannot open alt mtx file");

    let alt_reader = BufReader::new(alt_reader);
    let ref_reader = File::open(params.ref_mtx.to_string()).expect("cannot open ref mtx file");
    
    let ref_reader = BufReader::new(ref_reader);
    let mut used_loci: HashSet<usize> = HashSet::new();
    let mut line_number = 0;
    let mut total_loci = 0;
    let mut total_cells = 0;
    let mut all_loci: HashSet<usize> = HashSet::new();
    let mut locus_cell_counts: HashMap<usize, [u32; 2]> = HashMap::new();
    let mut locus_umi_counts: HashMap<usize, [u32; 2]> = HashMap::new();
    let mut locus_counts: HashMap<usize, HashMap<usize, [u32; 2]>> = HashMap::new();
    let mut cell_count_per_locus: HashMap<usize, usize> = HashMap::new(); // New hashmap Reza
    for (alt_line, ref_line) in izip!(alt_reader.lines(), ref_reader.lines()) {
        let alt_line = alt_line.expect("cannot read alt mtx");
        let ref_line = ref_line.expect("cannot read ref mtx");
        if line_number > 2 {
            let alt_tokens: Vec<&str> = alt_line.split_whitespace().collect();
            let ref_tokens: Vec<&str> = ref_line.split_whitespace().collect();
            let locus = alt_tokens[0].to_string().parse::<usize>().unwrap() - 1;
            all_loci.insert(locus);
            let cell = alt_tokens[1].to_string().parse::<usize>().unwrap() - 1;
            let ref_count = ref_tokens[2].to_string().parse::<u32>().unwrap();
            let alt_count = alt_tokens[2].to_string().parse::<u32>().unwrap();
            assert!(locus < total_loci);
            assert!(cell < total_cells);
            let cell_counts = locus_cell_counts.entry(locus).or_insert([0; 2]);
            let umi_counts = locus_umi_counts.entry(locus).or_insert([0; 2]);
            if ref_count > 0 { cell_counts[0] += 1; umi_counts[0] += ref_count; }
            if alt_count > 0 { cell_counts[1] += 1; umi_counts[1] += alt_count; }
            let cell_counts = locus_counts.entry(locus).or_insert(HashMap::new());
            cell_counts.insert(cell, [ref_count, alt_count]);
            *cell_count_per_locus.entry(locus).or_insert(0) += 1; // New hashmap Reza
        } else if line_number == 2 {
            let tokens: Vec<&str> = alt_line.split_whitespace().collect();
            total_loci = tokens[0].to_string().parse::<usize>().unwrap();
            total_cells = tokens[1].to_string().parse::<usize>().unwrap();
        }
        line_number += 1;
        
    }
    
    //REZA ref&alt counts per loci
    let mut ref_alt_counts_per_locus: HashMap<usize, [u32; 2]> = HashMap::new();
    for (&locus, cell_counts) in &locus_counts {
        let mut total_ref = 0;
        let mut total_alt = 0;
        for &counts in cell_counts.values() {
            total_ref += counts[0];
            total_alt += counts[1];
        }
        ref_alt_counts_per_locus.insert(locus, [total_ref, total_alt]);
    }


    let mut all_loci2: Vec<usize> = Vec::new();
    for loci in all_loci {
        all_loci2.push(loci);
    }
    let mut all_loci = all_loci2;

    all_loci.sort();
    let mut index_to_locus: Vec<usize> = Vec::new();
    let mut locus_to_index: HashMap<usize, usize> = HashMap::new();
    let mut cell_data: Vec<CellData> = Vec::new();
    for _cell in 0..total_cells {
        cell_data.push(CellData::new());
    }
    let mut locus_index = 0;
    for locus in all_loci {
        let cell_counts = locus_cell_counts.get(&locus).unwrap();
        let umi_counts = locus_umi_counts.get(&locus).unwrap();
        if cell_counts[0] >= params.min_ref && cell_counts[1] >= params.min_alt && umi_counts[0] >= params.min_ref_umis && umi_counts[1] >= params.min_alt_umis {
            used_loci.insert(locus);
            index_to_locus.push(locus);
            locus_to_index.insert(locus, locus_index);
            for (cell, counts) in locus_counts.get(&locus).unwrap() {
                if counts[0]+counts[1] == 0 { continue; }
                cell_data[*cell].alt_counts.push(counts[1]);
                cell_data[*cell].ref_counts.push(counts[0]);
                cell_data[*cell].loci.push(locus_index);
                cell_data[*cell].allele_fractions.push((counts[1] as f32)/((counts[0] + counts[1]) as f32));
                cell_data[*cell].log_binomial_coefficient.push(
                     statrs::function::factorial::ln_binomial((counts[1]+counts[0]) as u64, counts[1] as u64) as f32);
                cell_data[*cell].total_alleles += (counts[0] + counts[1]) as f32;
            
            }
            locus_index += 1;
        }
    }
    eprintln!("total loci used {}",used_loci.len());
    
    
    (used_loci.len(), total_cells, cell_data, index_to_locus, locus_to_index, used_loci, cell_count_per_locus, ref_alt_counts_per_locus)
}




struct CellData {
    allele_fractions: Vec<f32>,
    log_binomial_coefficient: Vec<f32>,
    alt_counts: Vec<u32>,
    ref_counts: Vec<u32>, // 
    loci: Vec<usize>,// ID
    total_alleles: f32,
}

impl CellData {
    fn new() -> CellData {
        CellData{
            allele_fractions: Vec::new(),
            log_binomial_coefficient: Vec::new(),
            alt_counts: Vec::new(),
            ref_counts: Vec::new(),
            loci: Vec::new(),
            total_alleles: 0.0,
        }
    }
}



fn load_barcodes(params: &Params) -> Vec<String> {

    let reader = reader(&params.barcodes);
    let mut cell_barcodes: Vec<String> = Vec::new();
    for line in reader.lines() {
        let line = line.expect("Unable to read line");
        cell_barcodes.push(line.to_string());
    }
    cell_barcodes
}


#[derive(Clone)]
struct Params {
    ref_mtx: String,
    alt_mtx: String,
    barcodes: String,
    num_clusters: usize,
    min_alt: u32,
    min_ref: u32,
    min_alt_umis: u32,
    min_ref_umis: u32,
    restarts: u32,
    known_cell_assignments: Option<String>,
    known_genotypes: Option<String>,
    known_genotypes_sample_names: Vec<String>,
    initialization_strategy: ClusterInit,
    threads: usize,
    seed: u8,
    cell_hashing_assignments: Option<String>,

}

#[derive(Clone)]
enum ClusterInit {
    KmeansPP,
    RandomUniform,
    RandomAssignment,
    MiddleVariance,
}

fn load_params() -> Params {
    let yaml = load_yaml!("params.yml");
    let params = App::from_yaml(yaml).get_matches();
    let ref_mtx = params.value_of("ref_matrix").unwrap();
    let alt_mtx = params.value_of("alt_matrix").unwrap();
    let barcodes = params.value_of("barcodes").unwrap();
    let num_clusters = params.value_of("num_clusters").unwrap();
    let num_clusters = num_clusters.to_string().parse::<usize>().unwrap();
    let min_alt = params.value_of("min_alt").unwrap_or("4");
    let min_alt = min_alt.to_string().parse::<u32>().unwrap();
    let min_ref = params.value_of("min_ref").unwrap_or("4");
    let min_ref = min_ref.to_string().parse::<u32>().unwrap();
    let restarts = params.value_of("restarts").unwrap_or("100");
    let restarts = restarts.to_string().parse::<u32>().unwrap();
    let known_cell_assignments = params.value_of("known_cell_assignments");
    let known_cell_assignments = match known_cell_assignments {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    let known_genotypes = params.value_of("known_genotypes");
    let known_genotypes = match known_genotypes {
        Some(x) => {
            assert!(known_cell_assignments == None, "Cannot set both known_genotypes and known_cell_assignments");
            Some(x.to_string())
        },
        None => None,
    };
    let known_genotypes_sample_names = params.values_of("known_genotypes_sample_names");
    let known_genotypes_sample_names: Vec<&str> = match known_genotypes_sample_names {
        Some(x) => x.collect(),
        None => Vec::new(),
    };
    let mut sample_names: Vec<String> = Vec::new();
    for name in known_genotypes_sample_names {
        sample_names.push(name.to_string());
    }

    let initialization_strategy = params.value_of("initialization_strategy").unwrap_or("random_uniform");
    let initialization_strategy = match initialization_strategy {
        "kmeans++" => ClusterInit::KmeansPP,
        "random_uniform" => ClusterInit::RandomUniform,
        "random_cell_assignment" => ClusterInit::RandomAssignment,
        "middle_variance" => ClusterInit::MiddleVariance,
        _ => {
            assert!(false, "initialization strategy must be one of kmeans++, random_uniform, random_cell_assignment, middle_variance");
            ClusterInit::RandomAssignment
        },
    };

    let threads = params.value_of("threads").unwrap_or("1");
    let threads = threads.to_string().parse::<usize>().unwrap();

    let seed = params.value_of("seed").unwrap_or("4");
    let seed = seed.to_string().parse::<u8>().unwrap();

    let min_ref_umis = params.value_of("min_ref_umis").unwrap_or("0");
    let min_ref_umis = min_ref_umis.to_string().parse::<u32>().unwrap();
    
    
    let min_alt_umis = params.value_of("min_alt_umis").unwrap_or("0");
    let min_alt_umis = min_alt_umis.to_string().parse::<u32>().unwrap();

    //reza
    let cell_hashing_assignments = params.value_of("cell_hashing_assignments");
    let cell_hashing_assignments = match cell_hashing_assignments {
        Some(x) => Some(x.to_string()),
        None => None,
    };
    //reza

    Params{
        ref_mtx: ref_mtx.to_string(),
        alt_mtx: alt_mtx.to_string(),
        barcodes: barcodes.to_string(),
        num_clusters: num_clusters,
        min_alt: min_alt,
        min_ref: min_ref,
        restarts: restarts,
        known_cell_assignments: known_cell_assignments,
        known_genotypes: known_genotypes,
        known_genotypes_sample_names: sample_names,
        initialization_strategy: initialization_strategy,
        threads: threads,
        seed: seed,
        min_alt_umis: min_alt_umis,
        min_ref_umis: min_ref_umis,
        cell_hashing_assignments: cell_hashing_assignments,
    }
}

fn new_seed(rng: &mut StdRng) -> [u8; 32] {
    let mut seed = [0; 32];
    for i in 0..32 {
        seed[i] = rng.gen::<u8>();
    }
    seed
}






fn total_loglikelihood_cell_hashing(
    cell_data: &Vec<CellData>,
    cluster_centers: &Vec<Vec<f32>>,
    log_prior: f32,
    cell_id_to_label_map: &HashMap<usize, String>)  -> f32 {

    let mut log_binom_loss = 0.0;
    let mut cell_log_loss: Vec<f32> = Vec::new();    

    for (celldex, cell) in cell_data.iter().enumerate() {

        let log_binoms = robust_binomial_loss(cell, &cluster_centers, log_prior, celldex);
        let log_loss = log_sum_exp(&log_binoms);
        log_binom_loss += log_loss;
        cell_log_loss.push(log_loss);
    }

    log_binom_loss 


}






fn robust_binomial_loss(cell_data: &CellData, cluster_centers: &Vec<Vec<f32>>, log_prior: f32, cellnum: usize) -> Vec<f32> {
    let mut log_probabilities: Vec<f32> = Vec::new();
    let mut sum = 0.0;
    
    for (cluster, center) in cluster_centers.iter().enumerate() {
         log_probabilities.push(log_prior);
        for (locus_index, locus) in cell_data.loci.iter().enumerate() {
            let alt_term = (cell_data.alt_counts[locus_index] as f32) * safe_ln(center[*locus]);
            let ref_term = (cell_data.ref_counts[locus_index] as f32) * safe_ln(1.0 - center[*locus]);
            log_probabilities[cluster] += cell_data.log_binomial_coefficient[locus_index] + alt_term + ref_term;
        }

        sum += log_probabilities[cluster];
    }
    
    log_probabilities
}

fn safe_ln(x: f32) -> f32 {
    if x > 0.0 {
        x.ln()
    } else {
        (1e-6f32).ln()
    }
}







//begining of reza snippet


fn print_cluster_centers_column_by_column(cluster_centers: &Vec<Vec<f32>>) {
    if cluster_centers.is_empty() {
        return;
    }

    let num_columns = cluster_centers[0].len();
    for col in 0..num_columns {
        println!("Column {}:", col);
        for (cluster_index, cluster) in cluster_centers.iter().enumerate() {
            if let Some(value) = cluster.get(col) {
                println!("  Cluster {}: {}", cluster_index, value);
            } else {
                println!("  Cluster {}: Out of bounds", cluster_index);
            }
        }
        println!(); 
    }
}



fn initialize_cluster_centers_cell_hashing(loci_used: &HashSet<usize>, 
    cell_data_from_sc: &Vec<CellData>, barcodes: Vec<String>, params: &Params) 
     -> (Vec<Vec<f32>>, HashMap<usize, String>, HashMap<String, (usize, String)>, Vec<(String, String)>) {

    let columns = loci_used.len(); 
    let rows = 2;
    let mut cluster_centers = vec![vec![0.0; (columns)]; rows];

    // HAYNES
    let mut barcode_to_cell_id_and_hash: HashMap<String, (usize, String)> = HashMap::new();
    let mut cell_id_to_barcode_and_hash: Vec<(String, String)> = Vec::new();

    let mut loci_vec: Vec<usize> = loci_used.iter().cloned().collect();
    loci_vec.sort_unstable();

    let assignments_file = match params.cell_hashing_assignments {
        Some(ref filename) => File::open(filename)
            .expect("Unable to open cell_hashing_assignments.tsv"),
        None => panic!("cell_hashing_assignments file not specified"),
    };
    
    let assignments: Vec<String> = BufReader::new(assignments_file)
        .lines()
        .collect::<Result<_, _>>()
        .expect("Unable to read lines from cell_hashing_assignments.tsv");

    


    let alt_file = File::open(&params.alt_mtx.to_string()).expect("Unable to open alt.mtx");
    let alt_lines: Vec<String> = BufReader::new(alt_file).lines()
        .skip(3) // Skip the first three lines
        .map(|line| line.expect("Unable to read line"))
        .collect();
        
    let ref_file = File::open(&params.ref_mtx).expect("Unable to open ref.mtx");
    let ref_lines: Vec<String> = BufReader::new(ref_file).lines()
        .skip(3) // Skip the first three lines
        .map(|line| line.expect("Unable to read line"))
        .collect();

    
    let mut assignments_map: HashMap<String, String> = HashMap::new();

    for (id,line) in assignments.iter().enumerate() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() == 2 {
            let barcode = parts[0].to_string();
            let assignment = parts[1].to_string();
            
            // Check if the assignment ends with "C1" or "A"
            if assignment.ends_with("C1") {
                assignments_map.insert(barcode.to_string(), "C1".to_string());
                barcode_to_cell_id_and_hash.insert(barcode.to_string(), (id, "C1".to_string()));
                cell_id_to_barcode_and_hash.push((barcode, "C1".to_string()));
            } else if assignment.ends_with("A") {
                assignments_map.insert(barcode.to_string(), "A".to_string());
                barcode_to_cell_id_and_hash.insert(barcode.to_string(), (id, "A".to_string()));
                cell_id_to_barcode_and_hash.push((barcode, "A".to_string()));
            }
        }
    }

    let mut cell_id_to_label_map: HashMap<usize, String> = HashMap::new();

    for (index, barcode) in barcodes.iter().enumerate() {
        let cell_id = index + 0;
        if let Some(label) = assignments_map.get(barcode) {
            // println!("barcode {} label {} cell_id {}", barcode, label, cell_id);
            cell_id_to_label_map.insert(cell_id, label.clone()); 
        }
    }

    let mut cell_to_locus_to_alt: HashMap<usize, HashMap<usize, usize>> = HashMap::new();
    let mut locus_to_cell_id: HashMap<usize, Vec<usize>> = HashMap::new();

    for alt_line in alt_lines.iter() {
        let tokens: Vec<&str> = alt_line.split_whitespace().collect();
        let locus = tokens[0].parse::<usize>().unwrap() - 1; 
        let cell = tokens[1].parse::<usize>().unwrap() - 1;    
        let count = tokens[2].parse::<usize>().unwrap(); 
    
        cell_to_locus_to_alt.entry(cell)
            .or_insert_with(HashMap::new)
            .insert(locus , count); 

        locus_to_cell_id.entry(locus)
        .or_insert_with(Vec::new)
        .push(cell); 

    }

    let mut cell_to_locus_to_ref: HashMap<usize, HashMap<usize, usize>> = HashMap::new();

    for ref_line in ref_lines.iter() {
        let tokens: Vec<&str> = ref_line.split_whitespace().collect();
        let locus = tokens[0].parse::<usize>().unwrap() - 1; 
        let cell = tokens[1].parse::<usize>().unwrap() - 1;    
        let count = tokens[2].parse::<usize>().unwrap(); 
    
        cell_to_locus_to_ref.entry(cell)
            .or_insert_with(HashMap::new)
            .insert(locus , count); 
    }

    generate_cell_hashing_cluster_centers(&cell_id_to_label_map, 
        &cell_to_locus_to_alt, &cell_to_locus_to_ref, 
        &locus_to_cell_id, &loci_vec, &mut cluster_centers);



    // for row in cluster_centers.iter_mut() {
    //     for value in row.iter_mut() {
    //         if *value == 0.0 {
    //             *value = 0.0001;
    //         }
    //     }
    // }


    (cluster_centers, cell_id_to_label_map, barcode_to_cell_id_and_hash, cell_id_to_barcode_and_hash)
}  



fn generate_cell_hashing_cluster_centers(cell_id_to_label_map: &HashMap<usize, String>, 
    cell_to_locus_to_alt: &HashMap<usize, HashMap<usize, usize>>, 
    cell_to_locus_to_ref: &HashMap<usize, HashMap<usize, usize>>, 
    locus_to_cell_id: &HashMap<usize, Vec<usize>>, loci_vec: &Vec<usize>, cluster_centers: &mut Vec<Vec<f32>>) {

    let mut A_num: usize = 0;
    let mut A_denom: usize = 0;
    let mut C1_denom: usize = 0;
    let mut C1_num: usize = 0;

    let mut current_locus_index = 0;
    let mut c1_no_counts = 0;
    let mut c2_no_counts = 0;
    for &temp_loci in loci_vec.iter() {

        let actual_loci = temp_loci + 0;
        C1_num = 0;
        C1_denom = 0;
        A_num = 0;
        A_denom = 0;


        match locus_to_cell_id.get(&actual_loci) {
            Some(cells) => {
                //wanna print the length of cells
                // println!("loci id = {}, number of cells {}", actual_loci , cells.len());

                for &cell in cells.iter() {
                
                    let mut label = String::new();
        
                    if let Some(l) = cell_id_to_label_map.get(&cell) {
                        label = l.clone();
                    }
                    
                    if label == "C1" {
                                
                        let mut alt_count = cell_to_locus_to_alt.get(&cell).and_then(|locus_map| locus_map.get(&actual_loci)).copied();
                        let mut alt_value: usize = 0;
                        if let Some(value) = alt_count {
                            C1_num += value;
                            alt_value = value;
                        }
                        let mut ref_count =cell_to_locus_to_ref.get(&cell).and_then(|locus_map| locus_map.get(&actual_loci)).copied();
                        if let Some(value) = ref_count {
                            C1_denom += value + alt_value;
                        }
        
                        
                    }
                    else if label == "A" {
                        let mut alt_count = cell_to_locus_to_alt.get(&cell).and_then(|locus_map| locus_map.get(&actual_loci)).copied();
                        
                        let mut alt_value: usize = 0;
                        if let Some(value) = alt_count {
                            A_num += value;
                            alt_value = value;
                        }
                        let mut ref_count =cell_to_locus_to_ref.get(&cell).and_then(|locus_map| locus_map.get(&actual_loci)).copied();
                        if let Some(value) = ref_count {
                            A_denom += value + alt_value;
                        }
                    }
            }

            },
            None => {
                println!("No cells found for locus ID {}", temp_loci);
            }

        }
        if C1_denom != 0 {
            cluster_centers[0][current_locus_index] = C1_num as f32 / C1_denom as f32
        } else {
            c1_no_counts += 1;
            //eprintln!("locus {} cluster 0 has no counts",current_locus_index);
        } 
        
        if A_denom != 0 {
            cluster_centers[1][current_locus_index] = A_num as f32 / A_denom as f32;
        } else {
            c2_no_counts += 1;
            //eprintln!("locus {} cluster 1 has no counts",current_locus_index);
        }
        cluster_centers[0][current_locus_index] = cluster_centers[0][current_locus_index].min(0.99).max(0.01);
        cluster_centers[1][current_locus_index] = cluster_centers[1][current_locus_index].min(0.99).max(0.01);
        

        current_locus_index += 1;

    }

        // eprintln!("cluster 1 had {} no_counts_loci and cluster 2 had {} no_counts_loci",c1_no_counts, c2_no_counts);

}



// end of reza snippet




// this annealing schedule is totally heuristic
                // this annealing schedule is based on the idea that the more data a cell has, 
                // the smoother we need to make the probability space to avoid near 1.0 and 0.0 probabilities
                // until the search is closer to converging