extern crate bio;
extern crate rust_htslib;

//use bio::io::fasta;

//use rust_htslib::bcf;
use rust_htslib::bcf::*;
use rust_htslib::bcf::Reader;
//use rust_htslib::prelude::*;

//use rust_htslib::bcf::header::Header;
//use rust_htslib::bcf::header::HeaderRecord;
use rust_htslib::bcf::record::*;

use bio::io::bed;



mod range;


use std::collections::HashMap;

#[derive(Eq, PartialEq, Hash)]
enum HaplotypeSide { Left, Right }


#[derive(Eq, PartialEq, Clone)]
struct Diff {
    pos: u32,
    reference: Vec<u8>,
    alternative: Vec<u8>,
}

#[derive(Eq, PartialEq, Hash)]
struct HaplotypeId {
    sample_id: usize,
    side: HaplotypeSide
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

fn load_diffs(reader: &mut IndexedReader, sample_count: usize) -> (HashMap<HaplotypeId, Vec<Diff>>, u32) {
    let mut xs : HashMap<HaplotypeId, Vec<Diff>> = HashMap::new();
    let mut parsed_number: u32 = 0;
    for r in reader.records() {
        match r {
            Ok(mut record) => {
                let position = record.pos();
                let alleles = record.alleles();
                let reference = alleles[0].to_vec();
                let alternative = alleles[1].to_vec();
                let number_of_alleles = alleles.len();
                let genotypes = record.genotypes().unwrap();
                parsed_number = parsed_number + 1;

                if number_of_alleles == 2 {
                    let diff = Diff { pos : position, reference : reference.clone(), alternative : alternative.clone() };
                    for sample_id in 0..sample_count {
                        let genotype = genotypes.get(sample_id);
                        assert!(number_of_alleles == genotype.len(), "Inconsistent number of alleles");

                        if has_alternative(&genotype, HaplotypeSide::Left) {
                            let haplotype_id = HaplotypeId { sample_id : sample_id, side: HaplotypeSide::Left };
                            xs.entry(haplotype_id).or_insert(vec![diff.clone()]).push(diff.clone());
                        }
                        if has_alternative(&genotype, HaplotypeSide::Right) {
                            let haplotype_id = HaplotypeId { sample_id : sample_id, side: HaplotypeSide::Right };
                            xs.entry(haplotype_id).or_insert(vec![diff.clone()]).push(diff.clone());
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

fn main() {
    let chromosome = "chr1";
    let bed_files: Vec<&str> = "Bcell-13,CD4-9,CD8-10,CLP-14,CMP-4,Erythro-15,GMP-5,HSC-1,LMPP-3,MCP,mDC,MEGA1,MEGA2,MEP-6,Mono-7,MPP-2,Nkcell-11,pDC".split(',').collect();
    let bcf = format!("/home/seb/masters/topmed/source/TOPMed_dbGaP_20180710/dbGaP-12336/65066/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.6a/phased/freeze.6a.{}.pass_only.phased.bcf", chromosome);

    let (merged_peaks, peak_map) = load_peak_files(&bed_files, chromosome);

    match IndexedReader::from_path(bcf) {
        Ok(mut reader) => {
            let rid = reader.header().name2rid(chromosome.as_bytes()).unwrap();
            let samples = reader.header().samples();
            let sample_count = samples.len();
            println!("Number of samples: {}", sample_count);
            for peak in merged_peaks {
                reader.fetch(rid, peak.start as u32, peak.end as u32).unwrap();
                let (xs, parsed_number) = load_diffs(&mut reader, sample_count);
                println!("Peak {} {} {}, {} variants", chromosome, peak.start, peak.end, parsed_number);
            }
        }
        Err(e) => println!("{}", e),
    }

}



#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}