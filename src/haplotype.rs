extern crate rust_htslib;

use rust_htslib::bcf::*;
use rust_htslib::bcf::record::*;

use super::types::*;
use super::util::*;
use super::range::*;
use std::collections::HashMap;
use std::rc::Rc;


fn load_diffs(reader: &mut IndexedReader, sample_positions_in_bcf: &Vec<usize>) -> (HashMap<HaplotypeId, Vec<Diff>>, u32) {
    let mut xs : HashMap<HaplotypeId, Vec<Diff>> = HashMap::new();
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
                    let mut sample_id = 0;
                    for &sample_position in sample_positions_in_bcf {
                        let genotype = genotypes.get(sample_position);
                        assert!(number_of_alleles == genotype.len(), "Inconsistent number of alleles");

                        let has_alternative_left = match unsafe { genotype.get_unchecked(0) } {
                            GenotypeAllele::Unphased(1) => true,
                            _ => false,
                        };
                        let has_alternative_right = match unsafe { genotype.get_unchecked(1) } {
                            GenotypeAllele::Phased(1) => true,
                            _ => false,
                        };
                        if has_alternative_left {
                            let haplotype_id = HaplotypeId { sample_id : sample_id, side: HaplotypeSide::Left };
                            xs.entry(haplotype_id).or_insert(Vec::new()).push(diff.clone());
                        }
                        if has_alternative_right {
                            let haplotype_id = HaplotypeId { sample_id : sample_id, side: HaplotypeSide::Right };
                            xs.entry(haplotype_id).or_insert(Vec::new()).push(diff.clone());
                        }
                        sample_id = sample_id + 1;
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

pub fn load_haplotypes(chromosome: &str, peak: &Range, reader: &mut IndexedReader, ref_haplotype: &Vec<NucleotidePos>, sample_positions_in_bcf: &Vec<usize>) -> (u32, HashMap<Vec<NucleotidePos>, Rc<Vec<HaplotypeId>>>) {
    let rid = reader.header().name2rid(chromosome.as_bytes()).unwrap();
    reader.fetch(rid, peak.start as u32, peak.end as u32).unwrap();
    let (xs, variant_count) = load_diffs(reader, sample_positions_in_bcf);
    let mut res = HashMap::new();
    for (diffs, haplotype_ids) in group_by_diffs(xs).drain(){
        let haplotype = patch_haplotype(peak, &diffs, ref_haplotype);
        res.insert(haplotype, haplotype_ids);
    }
    (variant_count, res)
}

fn get(r: Range, ref_genome_peak: &Vec<NucleotidePos>) -> Vec<NucleotidePos> {
    ref_genome_peak.iter().cloned().filter(|n| n.pos >= r.start && n.pos <= r.end).collect()
}

pub fn patch_haplotype(range: &Range, diffs: &Vec<Diff>, ref_haplotype: &Vec<NucleotidePos>) -> Vec<NucleotidePos> {
    let mut sorted_diffs: Vec<&Diff> = diffs.iter().filter(|d| d.pos >= range.start && d.pos <= range.end).collect();
    sorted_diffs.sort();

    fn next_chunk(range: &Range, ref_position: u64, ds: &Vec<&Diff>, ref_haplotype: &Vec<NucleotidePos>) -> Vec<NucleotidePos> {
        match ds.split_first() {
            None => {
                if ref_position > range.end {
                    return vec![];
                }
                else {
                    let chunk = get(Range::new(ref_position, range.end), ref_haplotype);
                    return chunk;
                }
            }
            Some((d, rest)) => {
                if d.pos > ref_position {
                    let mut chunk = get(Range::new(ref_position, d.pos-1), ref_haplotype);
                    chunk.extend(next_chunk(range, d.pos, ds, ref_haplotype).drain(..));
                    return chunk;
                }
                else if d.pos == ref_position && d.reference.len() == 1 { // SNV or insertion
                    let mut chunk = Vec::new();
                    //println!("{} {:#?} {:#?}", range, ref_haplotype, d);
                    //println!("assert {} == {}", d.reference[0], ref_haplotype[0].nuc);
                    let reference_first_nuc_at = {
                        let mut nuc = Nucleotide::N;
                        for np in ref_haplotype {
                            if np.pos == ref_position { nuc = np.nuc; }
                        }
                        nuc
                    };
                    if d.reference[0] != reference_first_nuc_at {
                        panic!("First reference nucleotide of variant doesn't match reference genome: {:#?}", d);
                    }
                    for i in &d.alternative {
                        chunk.push(NucleotidePos { pos: ref_position, nuc: i.clone() });
                    }
                    chunk.extend(next_chunk(range, ref_position + 1, &rest.to_vec(), ref_haplotype).drain(..));
                    return chunk;
                }
                else if d.pos == ref_position && d.alternative.len() == 1 { // Deletion
                    let mut chunk = vec![NucleotidePos { pos: ref_position, nuc: d.alternative[0].clone()}];
                    chunk.extend(next_chunk(range, ref_position + (d.reference.len() as u64), &rest.to_vec(), ref_haplotype).drain(..));
                    return chunk;
                }
                else if d.pos == ref_position {
                    panic!("Missing case in haplotype patcher");
                }
                else if ref_position >= range.end {
                    return get(Range::new(ref_position, ref_position), ref_haplotype);
                }
                else {
                    return vec![];
                }
            }

        }
    }

    return next_chunk(range, range.start, &sorted_diffs, ref_haplotype);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ref_haplotype() -> Vec<NucleotidePos> {
        vec![
            NucleotidePos { nuc: Nucleotide::A, pos: 0 },
            NucleotidePos { nuc: Nucleotide::C, pos: 1 },
            NucleotidePos { nuc: Nucleotide::G, pos: 2 },
            NucleotidePos { nuc: Nucleotide::T, pos: 3 }
        ]
    }


    #[test]
    fn test_patch_haplotype_with_no_diff() {
        let diffs = Vec::new();
        let patched = patch_haplotype(&Range::new(1,2), &diffs, &ref_haplotype());
        let expected = vec![nucp('C',1), nucp('G',2)];
        assert_eq!(patched, expected);

        let patched2 = patch_haplotype(&Range::new(0,2), &diffs, &ref_haplotype());
        let expected2 = vec![nucp('A',0), nucp('C',1), nucp('G',2)];
        assert_eq!(patched2, expected2);

        let patched3 = patch_haplotype(&Range::new(0,5), &diffs, &ref_haplotype());
        let expected3 = vec![nucp('A',0), nucp('C',1), nucp('G',2), nucp('T',3)];
        assert_eq!(patched3, expected3);
    }

    #[test]
    fn test_patch_haplotype_one_snp() {
        let diffs = vec![Diff { pos: 100, reference: nucs("A"), alternative: nucs("C") }];
        let patched = patch_haplotype(&Range::new(1,2), &diffs, &ref_haplotype());
        let expected = vec![nucp('C',1), nucp('G',2)];
        assert_eq!(patched, expected);

        let diffs2 = vec![Diff { pos: 1, reference: nucs("C"), alternative: nucs("N") }];
        let patched2 = patch_haplotype(&Range::new(1,2), &diffs2, &ref_haplotype());
        let expected2 = vec![nucp('N',1), nucp('G',2)];
        assert_eq!(patched2, expected2);

        let diffs3 = vec![Diff { pos: 2, reference: nucs("G"), alternative: nucs("A") }];
        let patched3 = patch_haplotype(&Range::new(1,2), &diffs3, &ref_haplotype());
        let expected3 = vec![nucp('C',1), nucp('A',2)];
        assert_eq!(patched3, expected3);
    }

    #[test]
    fn test_patch_haplotype_two_snp() {
        let diffs = vec![Diff { pos: 1, reference: nucs("C"), alternative: nucs("N")}, Diff { pos: 2, reference: nucs("G"), alternative: nucs("A") }];
        let patched = patch_haplotype(&Range::new(1,2), &diffs, &ref_haplotype());
        let expected = vec![nucp('N',1), nucp('A',2)];
        assert_eq!(patched, expected);

        let diffs2 = vec![Diff { pos: 1, reference: nucs("C"), alternative: nucs("N") }, Diff { pos: 4, reference: nucs("G"), alternative: nucs("A") }];
        let patched2 = patch_haplotype(&Range::new(1,2), &diffs2, &ref_haplotype());
        let expected2 = vec![nucp('N',1), nucp('G',2)];
        assert_eq!(patched2, expected2);
    }

    #[test]
    fn test_patch_haplotype_one_insert() {
        let diffs = vec![Diff { pos: 1, reference: nucs("C"), alternative: nucs("NN") }];
        let patched = patch_haplotype(&Range::new(1,2), &diffs, &ref_haplotype());
        let expected = vec![nucp('N',1), nucp('N',1), nucp('G',2)];
        assert_eq!(patched, expected);

        let diffs2 = vec![Diff { pos: 2, reference: nucs("G"), alternative: nucs("NN") }];
        let patched2 = patch_haplotype(&Range::new(1,2), &diffs2, &ref_haplotype());
        let expected2 = vec![nucp('C',1), nucp('N',2), nucp('N',2)];
        assert_eq!(patched2, expected2);

        let diffs3 = vec![Diff { pos: 3, reference: nucs("T"), alternative: nucs("NN") }];
        let patched3 = patch_haplotype(&Range::new(1,2), &diffs3, &ref_haplotype());
        let expected3 = vec![nucp('C',1), nucp('G',2)];
        assert_eq!(patched3, expected3);
    }

    #[test]
    fn test_patch_haplotype_one_deletion() {
        let diffs = vec![Diff { pos: 1, reference: nucs("CG"), alternative: nucs("C") }];
        let patched = patch_haplotype(&Range::new(1,2), &diffs, &ref_haplotype());
        let expected = vec![nucp('C',1)];
        assert_eq!(patched, expected);

        let diffs2 = vec![Diff { pos: 2, reference: nucs("GT"), alternative: nucs("G") }];
        let patched2 = patch_haplotype(&Range::new(1,2), &diffs2, &ref_haplotype());
        let expected2 = vec![nucp('C',1), nucp('G',2)];
        assert_eq!(patched2, expected2);

        // For simplicity, do not apply a Diff that starts before the window we're observing
        let diffs3 = vec![Diff { pos: 0, reference: nucs("AC"), alternative: nucs("A") }];
        let patched3 = patch_haplotype(&Range::new(1,2), &diffs3, &ref_haplotype());
        let expected3 = vec![nucp('C',1), nucp('G',2)];
        assert_eq!(patched3, expected3);
    }

    fn nucs(s: &str) -> Vec<Nucleotide> {
        let mut v = Vec::new();
        for c in s.chars() {
            if c == 'N' { v.push(Nucleotide::N); }
            else if c == 'A' { v.push(Nucleotide::A); }
            else if c == 'C' { v.push(Nucleotide::C); }
            else if c == 'G' { v.push(Nucleotide::G); }
            else if c == 'T' { v.push(Nucleotide::T); }
            else { panic!("Wrong nucleotide sequence in tests"); }
        }
        return v;
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