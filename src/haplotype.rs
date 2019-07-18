extern crate rust_htslib;
extern crate bio;

use rust_htslib::bcf::*;
use rust_htslib::bcf::record::*;

use super::types::*;
use super::util::*;
use super::range::*;
use std::collections::HashMap;
use std::rc::Rc;


fn load_diffs(reader: &mut IndexedReader) -> (HashMap<HaplotypeId, Vec<Diff>>, u32) {
    let mut xs : HashMap<HaplotypeId, Vec<Diff>> = HashMap::new();
    let sample_count = reader.header().sample_count() as usize;
    let mut variant_count: u32 = 0;
    for r in reader.records() {
        match r {
            Ok(mut record) => {
                let position = record.pos();
                let alleles = record.alleles();
                let reference = to_nucleotides(&alleles[0].to_vec());
                let alternative = to_nucleotides(&alleles[1].to_vec());
                let number_of_alleles = alleles.len();
                let genotypes = record.genotypes().unwrap();
                variant_count = variant_count + 1;

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
    (xs, variant_count)
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

pub fn load_haplotypes<F>(chromosome: &str, peak: &Range, reader: &mut IndexedReader, ref_genome: F) -> (u32, HashMap<Vec<NucleotidePos>, Rc<Vec<HaplotypeId>>>) where F: Fn(Range) -> Vec<NucleotidePos> {
    let rid = reader.header().name2rid(chromosome.as_bytes()).unwrap();
    reader.fetch(rid, peak.start as u32, peak.end as u32).unwrap();
    let (xs, variant_count) = load_diffs(reader);
    let mut res = HashMap::new();
    for (diffs, haplotype_ids) in group_by_diffs(xs).drain(){
        let haplotype = patch_haplotype(peak, &diffs, &ref_genome);
        res.insert(haplotype, haplotype_ids);
    }
    (variant_count, res)
}

pub fn patch_haplotype<F>(range: &Range, diffs: &Vec<Diff>, get: &F) -> Vec<NucleotidePos> where F: Fn(Range) -> Vec<NucleotidePos> {
    let mut sorted_diffs: Vec<&Diff> = diffs.iter().filter(|d| d.pos >= range.start && d.pos <= range.end).collect();
    sorted_diffs.sort();

    fn next_chunk<F>(range: &Range, ref_position: u64, ds: &Vec<&Diff>, get: &F) -> Vec<NucleotidePos> where F: Fn(Range) -> Vec<NucleotidePos> {
        match ds.split_first() {
            None => {
                if ref_position > range.end {
                    return vec![];
                }
                else {
                    let chunk = get(Range::new(ref_position, range.end));
                    return chunk;
                }
            }
            Some((d, rest)) => {
                if d.pos > ref_position {
                    let mut chunk = get(Range::new(ref_position, d.pos-1));
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
                    return get(Range::new(ref_position, ref_position));
                }
                else {
                    return vec![];
                }
            }

        }
    }

    return next_chunk(range, range.start, &sorted_diffs, get);
}