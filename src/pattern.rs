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

pub fn parse_pwm_files(pwm_file: &str, threshold_dir: &str, pwm_threshold: f32, wanted_pwms: Vec<String>, add_reverse_patterns: bool) -> Vec<Pattern> {

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
                                let pwm = Pattern::PWM { weights: weights.clone(), name: name.clone(), pattern_id: pattern_id, min_score: min_score, direction: PWMDirection::P};
                                pwms.push(pwm.clone());
                                if add_reverse_patterns {
                                    pwms.push(Pattern::PWM { weights: reverse_complement(weights.clone()), name: name.clone(), pattern_id: pattern_id, min_score: min_score, direction: PWMDirection::N });
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

fn apply_pwm(pattern: &Pattern, haplotype: &[NucleotidePos]) -> i32 {
    match pattern {
        Pattern::PWM{weights, name:_, pattern_id:_ , min_score:_ , direction:_} => {
            return weights.iter().zip(haplotype.iter()).map(|(weight, nucleotide)| select_weight(weight,nucleotide)).sum();
        }
        Pattern::OtherPattern{name:_, pattern_id:_} => {
            return 0;
        }
    }

}

//fn display_vec_nuc(vec: &Vec<Nucleotide>) -> String {
//    vec.iter().fold(String::new(), |acc, &arg| acc + &arg.to_string())
//}

pub fn matches(pattern: &Pattern, haplotype: &Vec<NucleotidePos>, haplotype_ids: Rc<Vec<HaplotypeId>>) -> Vec<Match> {
    match pattern {
        Pattern::PWM{weights, name:_, pattern_id, min_score, direction:_} => {
            let mut res = Vec::new();
            //println!("enter matches {} {}", haplotype.len(), pwm.weights.len());
            if haplotype.len() >= weights.len() {
                //println!("haplotype.len() >= pwm.weights.len() {} {} {}", haplotype.len(), pwm.weights.len(), pwm.name);
                for i in 0..(haplotype.len()-weights.len()+1) {
                    let score = apply_pwm(&pattern, &haplotype[i..]);
                    if score > *min_score {
                        //println!("score {} min_score {} name {} position {} direction {}", score, pwm.min_score, pwm.name, haplotype[i].pos, pwm.direction);
                        //println!("----- score {} min_score {}", score, pwm.min_score);
                        let m = Match { pos : haplotype[i].pos, pattern_id : *pattern_id, haplotype_ids: haplotype_ids.clone() };
                        res.push(m);
                    }
                }
            }
            else {
                //println!("Haplotype too short for pwm {} vs {}", haplotype.len(), pwm.weights.len());
            }
            return res;
        }
        Pattern::OtherPattern{name:_, pattern_id:_} => {
            return Vec::new();
        }
    }

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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_pwm_files() {
        let pwms = parse_pwm_files("HOCOMOCOv11_full_pwms_HUMAN_mono.txt", "thresholds", 0.001, vec!["GATA1_HUMAN.H11MO.1.A".to_string(), "GATA2_HUMAN.H11MO.1.A".to_string()], true);
        assert_eq!(pwms.len(), 4);
        assert_eq!(pwms[0], Pattern::PWM {weights: vec![Weight { acgtn: [322, -754, 193, -65, 0].to_vec() },
                                                        Weight { acgtn: [-490, 565, 200, -898, 0].to_vec() },
                                                        Weight { acgtn: [1022, -2694, -3126, 105, 0].to_vec() },
                                                        Weight { acgtn: [-4400, -4400, 1375, -3903, 0].to_vec() },
                                                        Weight { acgtn: [1377, -4400, -4400, -4400, 0].to_vec() },
                                                        Weight { acgtn: [-3325, -3126, -4400, 1363, 0].to_vec() },
                                                        Weight { acgtn: [1347, -3126, -3325, -2584, 0].to_vec() },
                                                        Weight { acgtn: [1296, -3573, -1421, -2584, 0].to_vec() },
                                                        Weight { acgtn: [-570, -357, 969, -2311, 0].to_vec() },
                                                        Weight { acgtn: [393, -220, 304, -1022, 0].to_vec() },
                                                        Weight { acgtn: [304, -144, 250, -705, 0].to_vec() }],
                                          name: "GATA1_HUMAN.H11MO.1.A".to_string(),
                                          pattern_id: 0,
                                          min_score: 4683,
                                          direction: PWMDirection::P});

        assert_eq!(pwms[1], Pattern::PWM {weights: vec![Weight { acgtn: [-705, 250, -144, 304, 0].to_vec() },
                                                        Weight { acgtn: [-1022, 304, -220, 393, 0].to_vec() },
                                                        Weight { acgtn: [-2311, 969, -357, -570, 0].to_vec() },
                                                        Weight { acgtn: [-2584, -1421, -3573, 1296, 0].to_vec() },
                                                        Weight { acgtn: [-2584, -3325, -3126, 1347, 0].to_vec() },
                                                        Weight { acgtn: [1363, -4400, -3126, -3325, 0].to_vec() },
                                                        Weight { acgtn: [-4400, -4400, -4400, 1377, 0].to_vec() },
                                                        Weight { acgtn: [-3903, 1375, -4400, -4400, 0].to_vec() },
                                                        Weight { acgtn: [105, -3126, -2694, 1022, 0].to_vec() },
                                                        Weight { acgtn: [-898, 200, 565, -490, 0].to_vec() },
                                                        Weight { acgtn: [-65, 193, -754, 322, 0].to_vec() }],
                                          name: "GATA1_HUMAN.H11MO.1.A".to_string(),
                                          pattern_id: 0,
                                          min_score: 4683,
                                          direction: PWMDirection::N});

        assert_eq!(pwms[2], Pattern::PWM {weights: vec![Weight { acgtn: [333, -754, 281, -210, 0].to_vec() },
                                                        Weight { acgtn: [-415, 551, 327, -1525, 0].to_vec() },
                                                        Weight { acgtn: [1093, -2961, -3325, -74, 0].to_vec() },
                                                        Weight { acgtn: [-4400, -3903, 1371, -3573, 0].to_vec() },
                                                        Weight { acgtn: [1355, -2694, -3325, -3903, 0].to_vec() },
                                                        Weight { acgtn: [-2584, -1770, -1600, 1268, 0].to_vec() },
                                                        Weight { acgtn: [1229, -1561, -2034, -1421, 0].to_vec() },
                                                        Weight { acgtn: [1117, -2311, -291, -2311, 0].to_vec() },
                                                        Weight { acgtn: [-516, -40, 814, -1681, 0].to_vec() },
                                                        Weight { acgtn: [509, -357, 388, -1818, 0].to_vec() },
                                                        Weight { acgtn: [509, -543, 91, -415, 0].to_vec() }],
                                          name: "GATA2_HUMAN.H11MO.1.A".to_string(),
                                          pattern_id: 1,
                                          min_score: 5314,
                                          direction: PWMDirection::P});


        assert_eq!(pwms[3], Pattern::PWM {weights: vec![Weight { acgtn: [-415, 91, -543, 509, 0].to_vec() },
                                                        Weight { acgtn: [-1818, 388, -357, 509, 0].to_vec() },
                                                        Weight { acgtn: [-1681, 814, -40, -516, 0].to_vec() },
                                                        Weight { acgtn: [-2311, -291, -2311, 1117, 0].to_vec() },
                                                        Weight { acgtn: [-1421, -2034, -1561, 1229, 0].to_vec() },
                                                        Weight { acgtn: [1268, -1600, -1770, -2584, 0].to_vec() },
                                                        Weight { acgtn: [-3903, -3325, -2694, 1355, 0].to_vec() },
                                                        Weight { acgtn: [-3573, 1371, -3903, -4400, 0].to_vec() },
                                                        Weight { acgtn: [-74, -3325, -2961, 1093, 0].to_vec() },
                                                        Weight { acgtn: [-1525, 327, 551, -415, 0].to_vec() },
                                                        Weight { acgtn: [-210, 281, -754, 333, 0].to_vec() }],
                                          name: "GATA2_HUMAN.H11MO.1.A".to_string(),
                                          pattern_id: 1,
                                          min_score: 5314,
                                          direction: PWMDirection::N});
    }

    #[test]
    fn test_parse_threshold() {
        let threshold = parse_threshold_file(&"thresholds/GATA1_HUMAN.H11MO.0.A.thr".to_string(), 1e-04);
        assert_eq!(threshold, Some(7751));
    }

    #[test]
    fn test_matches() {
        let c = Weight::new(0, 1000, 0, 0);
        let g = Weight::new(0, 0, 1000, 0);
        let pwm = Pattern::PWM {weights: vec![c,g], name: "pwm".to_string(), pattern_id: 5, min_score: 1500, direction: PWMDirection::P};
        let haplotype = vec![
            NucleotidePos { nuc: Nucleotide::A, pos: 10 },
            NucleotidePos { nuc: Nucleotide::C, pos: 11 },
            NucleotidePos { nuc: Nucleotide::G, pos: 12 },
            NucleotidePos { nuc: Nucleotide::T, pos: 13 }
        ];
        let haplotype_ids = Rc::new(Vec::new());
        let m = matches(&pwm, &haplotype, haplotype_ids.clone());
        let expected = vec![Match {pos: 11, pattern_id: 5, haplotype_ids: haplotype_ids.clone()}];
        assert_eq!(m, expected);
    }

    #[test]
    fn test_match_gataa() {
        let pwm = Pattern::PWM { weights: vec![Weight::new(0,0,100,0), Weight::new(100,0,0,0), Weight::new(0,0,0,100), Weight::new(100,0,0,0), Weight::new(100,0,0,0), ], name: "Example".to_string(), pattern_id: 123, min_score: 499, direction: PWMDirection::P };
        let haplotype_with_padding = vec![nucp('N',0), nucp('G',1), nucp('A',2), nucp('T',3), nucp('A',4), nucp('A',5), nucp('N',6)];
        let haplotype_without_padding = vec![nucp('G',1), nucp('A',2), nucp('T',3), nucp('A',4), nucp('A',5)];
        let haplotype_ids = Rc::new(vec![HaplotypeId { sample_id: 456, side: HaplotypeSide::Right }]);
        let ms = matches(&pwm, &haplotype_with_padding, haplotype_ids.clone());
        assert_eq!(ms.len(), 1);
        let ms2 = matches(&pwm, &haplotype_without_padding, haplotype_ids.clone());
        assert_eq!(ms2.len(), 1);

        let pwm2 = Pattern::PWM { weights: vec![Weight::new(0,0,100,0), Weight::new(100,0,0,0), Weight::new(0,0,0,100), Weight::new(100,0,0,0), Weight::new(100,0,0,0), ], name: "Example".to_string(), pattern_id: 123, min_score: 500, direction: PWMDirection::P };
        let ms3 = matches(&pwm2, &haplotype_with_padding, haplotype_ids.clone());
        assert_eq!(ms3.len(), 0);
        let ms4 = matches(&pwm2, &haplotype_without_padding, haplotype_ids.clone());
        assert_eq!(ms4.len(), 0);
    }

    fn nucp(c: char, p: u64) -> NucleotidePos {
        let nuc = {
            if c == 'N'      { Nucleotide::N }
            else if c == 'A' { Nucleotide::A }
            else if c == 'C' { Nucleotide::C }
            else if c == 'G' { Nucleotide::G }
            else if c == 'T' { Nucleotide::T }
            else { panic!("Wrong nucleotide sequence in tests"); }
        };
        return NucleotidePos { nuc: nuc, pos: p };
    }
}
