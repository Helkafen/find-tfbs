use std::collections::HashMap;
use std::io;
use std::fs::File;
use std::io::{BufRead};
use std::path::Path;
use std::rc::Rc;

use super::types::*;

pub fn parse_pwm_files(pwm_file: &str, threshold_file: &str, wanted_pwms: Vec<String>, add_reverse_patterns: bool) -> Vec<PWM> {
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
                                    let pwm = PWM { weights: current_weights.clone(), name: name.clone(), pattern_id: pattern_id, min_score: t, direction: PWMDirection::P };
                                    pwms.push(pwm);
                                    pattern_id = pattern_id + 1;
                                    {
                                        if add_reverse_patterns {
                                            let reverse_weights = { let mut x = current_weights; x.reverse(); x };
                                            let pwm = PWM { weights: reverse_weights, name: name, pattern_id: pattern_id, min_score: t, direction: PWMDirection::N };
                                            pwms.push(pwm);
                                            pattern_id = pattern_id + 1;
                                        }
                                    }
                                    current_weights = Vec::new();
                                },
                                None => println!("Couldn't find a PWM threshold for {}", name),
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
    pwms.into_iter().filter(|p| wanted_pwms.contains(&p.name)).collect()
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>> where P: AsRef<Path>, {
    let file = File::open(filename).expect("Could not open file");
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
    if haplotype.len() >= pwm.weights.len() {
        for i in 0..(haplotype.len()-pwm.weights.len()+1) {
            let score = apply_pwm(&pwm, &haplotype[i..]);
            if score > pwm.min_score {
                let m = Match { pos : haplotype[i].pos, pattern_id : pwm.pattern_id, haplotype_ids: haplotype_ids.clone() };
                res.push(m);
            }
        }
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
