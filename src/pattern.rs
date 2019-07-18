use std::collections::HashMap;
use std::io;
use std::fs::File;
use std::io::{BufRead};
use std::path::Path;
use std::rc::Rc;

use super::types::*;

pub fn parse_pwm_files(pwm_file: &str, threshold_file: &str) -> Vec<PWM> {
    fn parse_weight(s: &String) -> i32 {
        let x: f32 = s.parse().unwrap();
        (x * 1000.0).round() as i32
    }

    let mut thresholds = HashMap::new();

    if let Ok(lines) = read_lines(threshold_file) {
        for line in lines {
            if let Ok(l) = line {
                let x: Vec<String> = l.split_whitespace().into_iter().map(|a| a.to_string()).collect();
                if x.len() == 2 {
                    thresholds.insert(x[0].clone(), parse_weight(&x[1]));
                }
            }
        }
    }

    let mut pwms = Vec::new();
    let mut pattern_id = 0;
    let mut current_name: Option<String> = None;
    let mut current_weights: Vec<Weight> = Vec::new();
    if let Ok(lines) = read_lines(pwm_file) {
        for line in lines {
            if let Ok(l) = line {
                if l.chars().next() == Some('>') {
                    match current_name {
                        None => (),
                        Some(name) => {
                            match thresholds.get(&name) {
                                Some(&t) => {
                                    let pwm = PWM { weights: current_weights.into_iter().collect(), name: name, pattern_id: pattern_id, min_score: t };
                                    pwms.push(pwm);
                                    current_weights = Vec::new();
                                    pattern_id = pattern_id + 1;
                                },
                                None => println!("Couldn't find a PWM threshold"),
                            }

                        }
                    }
                    current_name = Some(l[1..].to_string());
                }
                else if l.is_empty() {
                }
                else {
                    let fields: Vec<String> = l.split_whitespace().into_iter().map(|a| a.to_string()).collect();
                    if fields.len() == 4 {
                        let w = Weight::new(parse_weight(&fields[0]), parse_weight(&fields[1]), parse_weight(&fields[2]), parse_weight(&fields[3]));
                        current_weights.push(w);
                    }
                }
            }
        }
    }
    match current_name {
        Some(name) => {
            match thresholds.get(&name) {
                Some(&t) => {
                    let pwm = PWM { weights: current_weights.into_iter().collect(), name: name, pattern_id: pattern_id, min_score: t };
                    pwms.push(pwm);
                }
                None => println!("Couldn't find a PWM threshold"),
            }
        }
        None => (),
    }

    pwms
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>> where P: AsRef<Path>, {
    let file = File::open(filename).expect("Could not open file"); // format!("Could not open file {}", filename)
    Ok(io::BufReader::new(file).lines())
}

fn apply_weight(w: &Weight, n: &NucleotidePos) -> i32 {
    match n.nuc {
        Nucleotide::A => w.w_a,
        Nucleotide::C => w.w_c,
        Nucleotide::G => w.w_g,
        Nucleotide::T => w.w_t,
        Nucleotide::N => 0,
    }
}

fn apply_weights(pwm: &PWM, haplotype: &[NucleotidePos]) -> i32 {
    return pwm.weights.iter().zip(haplotype.iter()).map(|(w, n)| apply_weight(w,n)).sum();
}

pub fn matches(pwm: &PWM, haplotype: &Vec<NucleotidePos>, haplotype_ids: Rc<Vec<HaplotypeId>>) -> Vec<Match> {
    let mut res = Vec::new();
    for i in 0..(haplotype.len()) {
        let score = apply_weights(&pwm, &haplotype[i..]);
        if score > pwm.min_score {
            let m = Match { pos : haplotype[i].pos, pattern_id : pwm.pattern_id, haplotype_ids: haplotype_ids.clone() };
            //println!("{} {} {} {}>{}",haplotype[i].pos, pwm.pattern_id, haplotype_ids.len(), score, pwm.min_score);
            res.push(m);
        }
    }
    return res;
}