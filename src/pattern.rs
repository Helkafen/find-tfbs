use std::collections::HashMap;
use std::io;
use std::fs::File;
use std::io::{BufRead};
use std::path::Path;
use std::rc::Rc;
use std::process::exit;
use std::fs::read_to_string;

use super::types::*;

fn parse_weight(s: &String) -> i32 {
    let x: f32 = s.parse().unwrap();
    (x * 1000.0).round() as i32
}

pub fn parse_threshold_file(filename: &String, pwm_threshold: f32) -> Option<i32> {
    let mut result = None;
    if let Ok(lines) = read_lines(filename.clone()) {
        for line in lines {
            if let Ok(l) = line {
                let x: Vec<String> = l.split_whitespace().into_iter().map(|a| a.to_string()).collect();
                if x.len() == 2 {
                    let weight = parse_weight(&x[0]);
                    let pvalue: f32 = x[1].parse().expect(&format!("Can't parse pvalue in file {}", &filename));
                    if pvalue > pwm_threshold {
                        result = Some(weight);
                    }
                }
            }
        }
    }
    return result;
}

pub fn parse_pwm_files(pwm_file: &str, threshold_dir: &str, pwm_threshold: f32, wanted_pwms: Vec<String>, add_reverse_patterns: bool) -> Vec<PWM> {

    let thresholds = {
        let mut thresholds = HashMap::new();
        for p in &wanted_pwms {
            let threshold_file = format!("{}/{}.thr", threshold_dir.trim_end_matches("/"), p);
            match parse_threshold_file(&threshold_file, pwm_threshold) {
                Some(min_score) => { thresholds.insert(p.clone(), min_score); },
                None => { println!("Could not parse {}",threshold_file); },
            };
        }
        thresholds
    };

    // See https://academic.oup.com/nar/article/46/D1/D252/4616875 and https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/nar/46/D1/10.1093_nar_gkx1106/2/gkx1106_supp.pdf for the definition of the PWMs
    // Careful: JDP2,GATA5,MESP1,ID4,ASCL2,SPIC have quality "D", which is too low to make de novo predictions
    // The pvalues are derived directly from the PWM matrix itself, following https://almob.biomedcentral.com/articles/10.1186/1748-7188-2-15
    // It is basically the proportion of haplotypes of the same length that match the PWM
    // At 0.0005, the median sensitivity to positive controls is 0.75
    // Maybe consider diPWMs, which seem to be a bit more specific: https://pdfs.semanticscholar.org/57e5/4745230282ff9692685dad6735bfe0d44942.pdf
    // Comparison of the ROC of several tools including ChipMunk: https://www.researchgate.net/publication/51632780_Tree-Based_Position_Weight_Matrix_Approach_to_Model_Transcription_Factor_Binding_Site_Profiles/figures?lo=1


    let mut pwms = Vec::new();
    let mut pattern_id = 0;

    match read_to_string(pwm_file) {
        Err(_) => { println!("Could not open file {}", pwm_file); exit(1); }
        Ok(content) => {
            for chunk in content.split(">") {
                if chunk.len() < 1 { continue; }
                let (name, weights) = parse_pwm_definition(chunk);
                    if wanted_pwms.contains(&name) {
                        match thresholds.get(&name) {
                            None => { println!("Couldn't find a PWM threshold for {}", name); },
                            Some(&min_score) => {
                                let pwm = PWM { weights: weights.clone(), name: name.clone(), pattern_id: pattern_id, min_score: min_score, direction: PWMDirection::P};
                                pwms.push(pwm.clone());
                                if add_reverse_patterns {
                                    pwms.push(PWM { weights: reverse_complement(weights.clone()), direction: PWMDirection::N, ..pwm });
                                }
                                println!("Loaded PWM {} (len {}, id {}, min_score {}) ", name, weights.len(), pattern_id, min_score);
                            }
                        }
                        pattern_id = pattern_id + 1;
                    }
                }
            }
        }
    pwms
}

fn parse_pwm_definition(chunk: &str) -> (String, Vec<Weight>) {
    let lines: Vec<&str> = chunk.split("\n").filter(|x| x.len() > 0).collect();
    let name = lines[0].to_string();
    let mut current_weights: Vec<Weight> = Vec::new();
    for line in lines.into_iter().skip(1) {
        let fields: Vec<String> = line.split_whitespace().into_iter().map(|a| a.to_string()).collect();
        if fields.len() == 4 {
            let w = Weight::new(parse_weight(&fields[0]), parse_weight(&fields[1]), parse_weight(&fields[2]), parse_weight(&fields[3]));
            current_weights.push(w);
        }
    }
    return (name, current_weights);
}

fn reverse_complement(w: Vec<Weight>) -> Vec<Weight> {
    let mut x = w;
    x.reverse();
    return x.iter().map(|w| complement(w)).collect();
}

fn complement(w: &Weight) -> Weight {
    let x = &w.acgtn;
    Weight::new(x[3], x[2], x[1], x[0])
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>> where P: AsRef<Path>, P: std::fmt::Display, P: std::clone::Clone {
    let file = File::open(filename.clone()).expect(&format!("Could not open file {}", filename));
    Ok(io::BufReader::new(file).lines())
}

fn select_weight(w: &Weight, n: &NucleotidePos) -> i32 {
    let i = n.nuc as usize;
    let x = unsafe { w.acgtn.get_unchecked(i) };
    *x
}

fn apply_pwm(pwm: &PWM, haplotype: &[NucleotidePos]) -> i32 {
    return pwm.weights.iter().zip(haplotype.iter()).map(|(weight, nucleotide)| select_weight(weight,nucleotide)).sum();
}

//fn display_vec_nuc(vec: &Vec<Nucleotide>) -> String {
//    vec.iter().fold(String::new(), |acc, &arg| acc + &arg.to_string())
//}

pub fn matches(pwm: &PWM, haplotype: &Vec<NucleotidePos>, haplotype_ids: Rc<Vec<HaplotypeId>>) -> Vec<Match> {
    let mut res = Vec::new();
    //println!("enter matches {} {}", haplotype.len(), pwm.weights.len());
    if haplotype.len() >= pwm.weights.len() {
        //println!("haplotype.len() >= pwm.weights.len() {} {} {}", haplotype.len(), pwm.weights.len(), pwm.name);
        for i in 0..(haplotype.len()-pwm.weights.len()+1) {
            let score = apply_pwm(&pwm, &haplotype[i..]);
            //println!("score {} min_score {} name {}", score, pwm.min_score, pwm.name);
            if score > pwm.min_score {
                //println!("----- score {} min_score {}", score, pwm.min_score);
                let m = Match { pos : haplotype[i].pos, pattern_id : pwm.pattern_id, haplotype_ids: haplotype_ids.clone() };
                res.push(m);
            }
        }
    }
    else {
        println!("Haplotype too short for pwm {} vs {}", haplotype.len(), pwm.weights.len());
    }
    return res;
}

            //println!("i:{} hap:{} pwm:{}", i, haplotype.len(), pwm.weights.len());
                //let mut matched = Vec::new();
                //for x in &haplotype[i..std::cmp::min(i+pwm.weights.len(), haplotype.len())] {
                //    matched.push(x.nuc);
                //}
                //println!("haplostart/end:{}/{} pos_of_match:{} {} {} {}>{} {}", haplotype[0].pos, haplotype[haplotype.len()-1].pos, haplotype[i].pos, pwm.pattern_id, haplotype_ids.len(), score, pwm.min_score, display_vec_nuc(&matched));
                //for w in &pwm.weights {
                //    print!("ACGTN ");
                //    for v in &w.acgtn {
                //        print!("{} ",v);
                //    }
                //    println!("");
                //}
