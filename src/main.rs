extern crate bio;
extern crate rust_htslib;
//extern crate itertools;

use bio::io::fasta;
use rust_htslib::bcf::*;
use rust_htslib::bcf::record::*;
use bio::io::bed;
use std::collections::HashMap;
use std::collections::HashSet;
use std::cmp;
use std::rc::Rc;
use std::path::Path;
use std::time::SystemTime;
//use itertools::Itertools;

mod range;



pub struct Weight {
    pub w_a : i32,
    pub w_c : i32,
    pub w_g : i32,
    pub w_t : i32,
    pub w_n : i32
}

impl Weight {
    pub fn new(a: i32, c: i32, g: i32, t: i32, n: i32) -> Weight {
        return Weight { w_a : a, w_c : c, w_g : g, w_t : t, w_n : n };
    }
}

pub struct PWM { pub weights: Vec<Weight>, pub name: String, pub pattern_id: u16, min_score: i32 }

#[derive(Eq, PartialEq, Debug)]
pub struct Match {
    pub pos: u64,
    pub pattern_id: u16,
    haplotype_ids: Rc<Vec<HaplotypeId>>
}

#[derive(Eq, PartialEq, Clone, Ord, PartialOrd, Debug, Hash)]
pub enum Nucleotide {
    A, C, G, T, N
}

#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Clone, Hash)]
pub struct NucleotidePos {
    pub nuc: Nucleotide,
    pub pos: u64
}

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
enum HaplotypeSide { Left, Right }

#[derive(Eq, PartialEq, Clone, Ord, PartialOrd, Hash, Debug)]
struct Diff {
    pos: u64,
    reference: Vec<Nucleotide>,
    alternative: Vec<Nucleotide>,
}

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub struct HaplotypeId {
    sample_id: usize,
    side: HaplotypeSide
}


fn apply_weight(w: &Weight, n: &NucleotidePos) -> i32 {
    match n.nuc {
        Nucleotide::A => w.w_a,
        Nucleotide::C => w.w_c,
        Nucleotide::G => w.w_g,
        Nucleotide::T => w.w_t,
        Nucleotide::N => w.w_n,
    }
}

pub fn apply_weights(pwm: &PWM, haplotype: &[NucleotidePos]) -> i32 {
    return pwm.weights.iter().zip(haplotype.iter()).map(|(w, n)| apply_weight(w,n)).sum();
}

pub fn matches(pwm: &PWM, haplotype: &Vec<NucleotidePos>, haplotype_ids: Rc<Vec<HaplotypeId>>) -> Vec<Match> {
    let mut res = Vec::new();
    for i in 0..(haplotype.len()) {
        let score = apply_weights(&pwm, &haplotype[i..]);
        if score > pwm.min_score {
            res.push(Match { pos : haplotype[i].pos, pattern_id : pwm.pattern_id, haplotype_ids: haplotype_ids.clone() });
        }
    }
    return res;
}

fn has_alternative(genotype: &Genotype, side: HaplotypeSide) -> bool{
    let index = match side {
        HaplotypeSide::Left => 0,
        HaplotypeSide::Right => 1,
    };
    match unsafe { genotype.get_unchecked(index) } {
        GenotypeAllele::Phased(1) => true,
        GenotypeAllele::Unphased(1) => true,
        _ => false,
    }
}

fn load_bed(filename: &String, chromosome: &str) -> Vec<range::Range> {
    let mut reader = bed::Reader::from_file(filename).unwrap();
    let mut xs = Vec::new();
    for record in reader.records() {
        let rec = record.unwrap();
        if rec.chrom() == chromosome {
            xs.push(range::Range::new(rec.start(), rec.end()));
        }
    }
    xs
}

fn sum_peak_sizes(peaks: &Vec<range::Range>) -> u64 {
    peaks.iter().map(|r| r.end-r.start).sum()
}

fn load_peak_files(bed_files: &Vec<&str>, chromosome: &str) -> (Vec<range::Range>, HashMap<String, Vec<range::Range>>) {
    let mut peak_map: HashMap<String, Vec<range::Range>> = HashMap::new();
    for bed_file in bed_files {
        let peaks = load_bed(&format!("/home/seb/masters/regu/dnamotifs/bed/{}.bed", bed_file), chromosome);
        println!("Loaded {}.bed:\t {} peaks covering {} bp", bed_file, peaks.len(), sum_peak_sizes(&peaks));
        peak_map.insert(bed_file.to_string(), peaks);
    }
    let vals: Vec<Vec<range::Range>> = peak_map.values().map(|x| x.clone()).collect();
    let rs: range::RangeStack = vals.concat().iter().collect();
    println!("Merged all peak files: {} peaks covering {} bp", rs.ranges.len(), sum_peak_sizes(&rs.ranges));
    (rs.ranges, peak_map)
}

fn to_nucleotide(l: u8) -> Nucleotide {
    if      l == 65 { return Nucleotide::A; }
    else if l == 67 { return Nucleotide::C; }
    else if l == 71 { return Nucleotide::G; }
    else if l == 84 { return Nucleotide::T; }
    else if l == 78 { return Nucleotide::N; }
    else if l == 97 { return Nucleotide::A; }
    else if l == 99 { return Nucleotide::C; }
    else if l == 103 { return Nucleotide::G; }
    else if l == 116 { return Nucleotide::T; }
    else if l == 110 { return Nucleotide::N; }
    else { panic!("Unknown nucleotide {}", l); }
}

fn to_nucleotides(letters: &Vec<u8>) -> Vec<Nucleotide> {
    return letters.iter().map(|&l| to_nucleotide(l)).collect();
}

fn to_nucleotides_pos(letters: &Vec<u8>, r: &range::Range) -> Vec<NucleotidePos> {
    let nucleotides = to_nucleotides(letters);
    let mut res = Vec::with_capacity(nucleotides.len());
    let mut pos = r.start;
    for n in nucleotides {
        res.push(NucleotidePos {nuc: n, pos: pos });
        pos = pos + 1;
    }
    res
}

fn patch_haplotype<F>(range: range::Range, diffs: &Vec<Diff>, get: F) -> Vec<NucleotidePos> where F: Fn(range::Range) -> Vec<NucleotidePos> {
    let mut sorted_diffs: Vec<&Diff> = diffs.iter().filter(|d| d.pos >= range.start && d.pos <= range.end).collect();
    sorted_diffs.sort();

    fn next_chunk<F>(range: range::Range, ref_position: u64, ds: &Vec<&Diff>, get: F) -> Vec<NucleotidePos> where F: Fn(range::Range) -> Vec<NucleotidePos> {
        match ds.split_first() {
            None => {
                if ref_position > range.end {
                    return vec![];
                }
                else {
                    let chunk = get(range::Range::new(ref_position, range.end));
                    return chunk;
                }
            }
            Some((d, rest)) => {
                if d.pos > ref_position {
                    let mut chunk = get(range::Range::new(ref_position, d.pos-1));
                    chunk.extend(next_chunk(range, d.pos, ds, get).drain(..));
                    return chunk;
                }
                else if d.pos == ref_position && d.reference.len() == 1 { // SNV or insertion
                    let mut chunk = Vec::new();
                    for i in &d.alternative {
                        chunk.push(NucleotidePos { pos: ref_position, nuc: i.clone() });
                    }
                    chunk.extend(next_chunk(range, ref_position + 1, &rest.to_vec(), get).drain(..));
                    return chunk;
                }
                else if d.pos == ref_position && d.alternative.len() == 1 { // Deletion
                    let mut chunk = vec![NucleotidePos { pos: ref_position, nuc: d.alternative[0].clone()}];
                    chunk.extend(next_chunk(range, ref_position + (d.reference.len() as u64), &rest.to_vec(), get).drain(..));
                    return chunk;
                }
                else if d.pos == ref_position {
                    panic!("Missing case in haplotype patcher");
                }
                else if ref_position >= range.end {
                    return get(range::Range::new(ref_position, ref_position));
                }
                else {
                    return vec![];
                }
            }

        }
    }

    return next_chunk(range, range.start, &sorted_diffs, get);
}


fn load_diffs(reader: &mut IndexedReader, sample_count: usize) -> (HashMap<HaplotypeId, Vec<Diff>>, u32) {
    let mut xs : HashMap<HaplotypeId, Vec<Diff>> = HashMap::new();
    let mut parsed_number: u32 = 0;
    for r in reader.records() {
        match r {
            Ok(mut record) => {
                let position = record.pos();
                let alleles = record.alleles();
                let reference = to_nucleotides(&alleles[0].to_vec());
                let alternative = to_nucleotides(&alleles[1].to_vec());
                let number_of_alleles = alleles.len();
                let genotypes = record.genotypes().unwrap();
                parsed_number = parsed_number + 1;

                if number_of_alleles == 2 {
                    let diff = Diff { pos : position as u64, reference : reference.clone(), alternative : alternative.clone() };
                    for sample_id in 0..sample_count {
                        let genotype = genotypes.get(sample_id);
                        assert!(number_of_alleles == genotype.len(), "Inconsistent number of alleles");

                        if has_alternative(&genotype, HaplotypeSide::Left) {
                            let haplotype_id = HaplotypeId { sample_id : sample_id, side: HaplotypeSide::Left };
                            xs.entry(haplotype_id).or_insert(Vec::new()).push(diff.clone());
                        }
                        if has_alternative(&genotype, HaplotypeSide::Right) {
                            let haplotype_id = HaplotypeId { sample_id : sample_id, side: HaplotypeSide::Right };
                            xs.entry(haplotype_id).or_insert(Vec::new()).push(diff.clone());
                        }
                    }
                }
                else {
                    println!("Unusual number of alleles: {}", number_of_alleles)
                }

            },
            Err(e) => println!("Bad record {}", e),
        }
    }
    (xs, parsed_number)
}


fn group_by_diffs(mut diffs: HashMap<HaplotypeId, Vec<Diff>>) -> HashMap<Vec<Diff>, Rc<Vec<HaplotypeId>>> {
    let mut res: HashMap<Vec<Diff>, Vec<HaplotypeId>> = HashMap::new();
    for (h,d) in diffs.drain() {
        res.entry(d).or_insert(Vec::new()).push(h);
    }
    let mut rc: HashMap<Vec<Diff>, Rc<Vec<HaplotypeId>>> = HashMap::new();
    for (k,v) in res.drain() {
        rc.insert(k, Rc::new(v));
    }
    return rc;
}

fn run_matches(i: u64) {
    let c = Weight::new(0, 1000, 0, 0, 0);
    let g = Weight::new(0, 0, 1000, 0, 0);
    let pwm = PWM {weights: vec![c,g], name: "pwm".to_string(), pattern_id: 5, min_score: 1500};
    let haplotype = vec![
        NucleotidePos { nuc: Nucleotide::A, pos: 10 },
        NucleotidePos { nuc: Nucleotide::C, pos: 11 },
        NucleotidePos { nuc: Nucleotide::G, pos: 12 },
        NucleotidePos { nuc: Nucleotide::T, pos: i }
    ];
    let haplotype_ids = Rc::new(Vec::new());
    let m = matches(&pwm, &haplotype, haplotype_ids.clone());
    let expected = vec![Match {pos: 11, pattern_id: 5, haplotype_ids: haplotype_ids.clone()}];
    if( i + m.len() as u64) % 10000000 == 0 {
        println!("{} {}",i, m.len());
    }
    assert_eq!(m, expected);
}

fn get_sample_names(reader: &IndexedReader) -> Vec<String> {
    let mut x = Vec::new();
    for i in reader.header().samples().iter() {
        x.push(String::from_utf8_lossy(i).into_owned());
    }
    x
}

fn repeat(size: usize, val: u32) -> Vec<u32> {
    let mut x = Vec::new();
    for _i in 0..size {
        x.push(val);
    }
    x
}

fn select_inner_peaks(peak: range::Range, peak_map: &HashMap<String, Vec<range::Range>>) -> HashMap<&String, Vec<&range::Range>> {
    let mut ip = HashMap::new();
    for (s,ps) in peak_map.iter() {
        for p in ps {
            if p.overlaps(&peak) {
                ip.entry(s).or_insert(Vec::new()).push(p)
            }
        }
    }
    ip
}


fn main() {
    let chromosome = "chr1";
    let bed_files: Vec<&str> = "Bcell-13,CD4-9,CD8-10,CLP-14,CMP-4,Erythro-15,GMP-5,HSC-1,LMPP-3,MCP,mDC,MEGA1,MEGA2,MEP-6,Mono-7,MPP-2,Nkcell-11,pDC".split(',').collect();
    let bcf = format!("/home/seb/masters/topmed/source/TOPMed_dbGaP_20180710/dbGaP-12336/65066/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.6a/phased/freeze.6a.{}.pass_only.phased.bcf", chromosome);
    let reference_genome_file = "/home/seb/masters/hg38.fa";

    let pwm_list: Vec<PWM> = Vec::new();

    let (merged_peaks, peak_map) = load_peak_files(&bed_files, chromosome);

    let mut reader = IndexedReader::from_path(bcf).expect("Error while opening the bcf file");
    let mut reference_genome = bio::io::fasta::IndexedReader::from_file(&Path::new(reference_genome_file)).expect("Error while opening the reference genome");

    let rid = reader.header().name2rid(chromosome.as_bytes()).unwrap();
    let samples = get_sample_names(&mut reader);
    let sample_count = samples.len();
    let null_count: Vec<u32> = repeat(sample_count, 0);
    println!("Number of samples: {}", sample_count);
    let mut peak_id = 0;
    let number_of_peaks = &merged_peaks.len();
    let all_haplotypes_with_reference_genome: HashSet<HaplotypeId> = {
            let mut x = HashSet::new();
            for i in (0..sample_count).into_iter() {
                x.insert(HaplotypeId {sample_id: i, side: HaplotypeSide::Left});
                x.insert(HaplotypeId {sample_id: i, side: HaplotypeSide::Right});
            }
            x
        };
    let start_time = SystemTime::now();

    for peak in merged_peaks {
        let peak_start_time = SystemTime::now();
        peak_id = peak_id + 1;
        reference_genome.fetch(chromosome, peak.start, peak.end).expect("Error while seeking in reference genome file");
        let ref_genome_peak: Vec<NucleotidePos> = {
            let mut text = Vec::new();
            reference_genome.read(&mut text).expect("Error while reading in reference genome file");
            to_nucleotides_pos(&text, &peak)
        };
        let ref_genome = |r: range::Range| -> Vec<NucleotidePos> {
            ref_genome_peak.iter().cloned().filter(|n| n.pos >= r.start && n.pos <= r.end).collect()
        };

        let inner_peaks: HashMap<&String, Vec<&range::Range>> = select_inner_peaks(peak, &peak_map);
        reader.fetch(rid, peak.start as u32, peak.end as u32).unwrap();
        let (xs, parsed_number) = load_diffs(&mut reader, sample_count);
        let mut match_list = Vec::new();
        let mut number_of_haplotypes = 0;
        let mut haplotypes_with_reference_genome: HashSet<HaplotypeId> = all_haplotypes_with_reference_genome.iter().cloned().collect();

        for (diffs, haplotype_ids) in group_by_diffs(xs).drain() {
            for h in haplotype_ids.iter()  {
                haplotypes_with_reference_genome.remove(&h);
            }
            let patched_haplotype = patch_haplotype(peak, &diffs, ref_genome);
            for pwm in &pwm_list {
                match_list.extend(matches(pwm, &patched_haplotype, haplotype_ids.clone()));
            }
            number_of_haplotypes = number_of_haplotypes + 1;
        };
        if !haplotypes_with_reference_genome.is_empty() {
            number_of_haplotypes = number_of_haplotypes + 1;
            let haplotype = ref_genome(peak.clone());
            let x: Vec<HaplotypeId> = haplotypes_with_reference_genome.into_iter().collect();
            let hap_ids = Rc::new(x);
            for pwm in &pwm_list {
                match_list.extend(matches(pwm, &haplotype, hap_ids.clone()));
            }
        }

        let mut counts = count_matches_by_sample(&match_list, &inner_peaks, &null_count);
        let fake_genotypes = {
            let mut x = HashMap::new();
            for (k,v) in counts.drain() {
                x.insert(k, counts_as_genotypes(v));
            }
            x
        };
        let number_of_matches: u64 = match_list.iter().map(|m| m.haplotype_ids.len() as u64).sum();
        let peak_time_elapsed = peak_start_time.elapsed().unwrap().as_millis();
        let global_time_elapsed = start_time.elapsed().unwrap().as_millis();
        println!("Peak {}/{}\t{} ms ({} total)\t{}\t{}\t{} haplotypes\t{} variants\t{} matches", peak_id, number_of_peaks, peak_time_elapsed, global_time_elapsed, peak.start, peak.end, number_of_haplotypes, parsed_number, number_of_matches);
    }

}

fn counts_as_genotypes(v: Vec<u32>) -> String {
    let mut res = String::with_capacity(v.len()*4);
    let min = v.iter().min();
    let max = v.iter().max();
    match (min, max) {
        (Some(&lowest), Some(&highest)) => {
            let intermediate = (lowest + highest) / 2;
            for x in v {
                if x == lowest {res.push_str("\t0|0");}
                else if x == highest {res.push_str("\t1|1");}
                else if x == intermediate {res.push_str("\t0|1");}
                else if x - lowest < intermediate - x {res.push_str("\t0|0");}
                else if highest - x < x - intermediate {res.push_str("\t1|1");}
                else {res.push_str("\t0|1");}
            }
        },
        (None, _) => (),
        (_, None) => (),
    }
    res
}

fn count_matches_by_sample<'a>(match_list: &Vec<Match>, inner_peaks: &'a HashMap<&String, Vec<&range::Range>>, null_count: &Vec<u32>) -> HashMap<(&'a String, &'a range::Range, u16), Vec<u32>> {
    let mut ppp = HashMap::new();
    for m in match_list {
        let pos = m.pos;
        for (&s,inner) in inner_peaks.iter().map(|(s,x)| (s, x.iter().filter(|y| y.contains(pos)))) {
            for &inner_peak in inner {
                let key = (s, inner_peak, m.pattern_id);
                for haplotype_id in m.haplotype_ids.iter() {
                    match ppp.entry(key).or_insert(null_count.clone()).get_mut(haplotype_id.sample_id as usize) {
                        None => panic!("Wrong sample_id index"),
                        Some(x) => *x = *x+1,
                    }
                }
            }
        }
    }
    ppp
}

fn ref_genome(r: range::Range) -> Vec<NucleotidePos> {
    let chromosome = vec![
        NucleotidePos { nuc: Nucleotide::A, pos: 0 },
        NucleotidePos { nuc: Nucleotide::C, pos: 1 },
        NucleotidePos { nuc: Nucleotide::G, pos: 2 },
        NucleotidePos { nuc: Nucleotide::T, pos: 3 }
        ];
    return chromosome[r.start as usize .. cmp::min((r.end+1) as usize,chromosome.len())].to_vec();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn test_matches() {
        let c = Weight::new(0, 1000, 0, 0, 0);
        let g = Weight::new(0, 0, 1000, 0, 0);
        let pwm = PWM {weights: vec![c,g], name: "pwm".to_string(), pattern_id: 5, min_score: 1500};
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
    fn test_patch_haplotype_with_no_diff() {
        let diffs = Vec::new();
        let patched = patch_haplotype(range::Range::new(1,2), &diffs, ref_genome);
        let expected = vec![NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched, expected);

        let patched2 = patch_haplotype(range::Range::new(0,2), &diffs, ref_genome);
        let expected2 = vec![NucleotidePos { nuc: Nucleotide::A, pos: 0 }, NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched2, expected2);

        let patched3 = patch_haplotype(range::Range::new(0,5), &diffs, ref_genome);
        let expected3 = vec![NucleotidePos { nuc: Nucleotide::A, pos: 0 }, NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }, NucleotidePos { nuc: Nucleotide::T, pos: 3 }];
        assert_eq!(patched3, expected3);
    }

    #[test]
    fn test_patch_haplotype_one_snp() {
        let diffs = vec![Diff { pos: 100, reference: vec![Nucleotide::A], alternative: vec![Nucleotide::C] }];
        let patched = patch_haplotype(range::Range::new(1,2), &diffs, ref_genome);
        let expected = vec![NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched, expected);

        let diffs2 = vec![Diff { pos: 1, reference: vec![Nucleotide::C], alternative: vec![Nucleotide::N] }];
        let patched2 = patch_haplotype(range::Range::new(1,2), &diffs2, ref_genome);
        let expected2 = vec![NucleotidePos { nuc: Nucleotide::N, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched2, expected2);

        let diffs3 = vec![Diff { pos: 2, reference: vec![Nucleotide::G], alternative: vec![Nucleotide::A] }];
        let patched3 = patch_haplotype(range::Range::new(1,2), &diffs3, ref_genome);
        let expected3 = vec![NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::A, pos: 2 }];
        assert_eq!(patched3, expected3);
    }

    #[test]
    fn test_patch_haplotype_two_snp() {
        let diffs = vec![Diff { pos: 1, reference: vec![Nucleotide::C], alternative: vec![Nucleotide::N] }, Diff { pos: 2, reference: vec![Nucleotide::G], alternative: vec![Nucleotide::A] }];
        let patched = patch_haplotype(range::Range::new(1,2), &diffs, ref_genome);
        let expected = vec![NucleotidePos { nuc: Nucleotide::N, pos: 1 }, NucleotidePos { nuc: Nucleotide::A, pos: 2 }];
        assert_eq!(patched, expected);

        let diffs2 = vec![Diff { pos: 1, reference: vec![Nucleotide::C], alternative: vec![Nucleotide::N] }, Diff { pos: 4, reference: vec![Nucleotide::G], alternative: vec![Nucleotide::A] }];
        let patched2 = patch_haplotype(range::Range::new(1,2), &diffs2, ref_genome);
        let expected2 = vec![NucleotidePos { nuc: Nucleotide::N, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched2, expected2);
    }

    #[test]
    fn test_patch_haplotype_one_insert() {
        let diffs = vec![Diff { pos: 1, reference: vec![Nucleotide::C], alternative: vec![Nucleotide::N, Nucleotide::N] }];
        let patched = patch_haplotype(range::Range::new(1,2), &diffs, ref_genome);
        let expected = vec![NucleotidePos { nuc: Nucleotide::N, pos: 1 }, NucleotidePos { nuc: Nucleotide::N, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched, expected);

        let diffs2 = vec![Diff { pos: 2, reference: vec![Nucleotide::C], alternative: vec![Nucleotide::N, Nucleotide::N] }];
        let patched2 = patch_haplotype(range::Range::new(1,2), &diffs2, ref_genome);
        let expected2 = vec![NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::N, pos: 2 }, NucleotidePos { nuc: Nucleotide::N, pos: 2 }];
        assert_eq!(patched2, expected2);

        let diffs3 = vec![Diff { pos: 3, reference: vec![Nucleotide::C], alternative: vec![Nucleotide::N, Nucleotide::N] }];
        let patched3 = patch_haplotype(range::Range::new(1,2), &diffs3, ref_genome);
        let expected3 = vec![NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched3, expected3);
    }

    #[test]
    fn test_patch_haplotype_one_deletion() {
        let diffs = vec![Diff { pos: 1, reference: vec![Nucleotide::C, Nucleotide::G], alternative: vec![Nucleotide::C] }];
        let patched = patch_haplotype(range::Range::new(1,2), &diffs, ref_genome);
        let expected = vec![NucleotidePos { nuc: Nucleotide::C, pos: 1 }];
        assert_eq!(patched, expected);

        let diffs2 = vec![Diff { pos: 2, reference: vec![Nucleotide::G, Nucleotide::T], alternative: vec![Nucleotide::G] }];
        let patched2 = patch_haplotype(range::Range::new(1,2), &diffs2, ref_genome);
        let expected2 = vec![NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched2, expected2);

        // For simplicity, do not apply a Diff that starts before the window we're observing
        let diffs3 = vec![Diff { pos: 0, reference: vec![Nucleotide::A, Nucleotide::C], alternative: vec![Nucleotide::A] }];
        let patched3 = patch_haplotype(range::Range::new(1,2), &diffs3, ref_genome);
        let expected3 = vec![NucleotidePos { nuc: Nucleotide::C, pos: 1 }, NucleotidePos { nuc: Nucleotide::G, pos: 2 }];
        assert_eq!(patched3, expected3);
    }
}