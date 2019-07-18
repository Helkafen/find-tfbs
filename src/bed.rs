use super::range;

use bio::io::bed;
use std::collections::HashMap;



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

pub fn load_peak_files(bed_files: &Vec<&str>, chromosome: &str) -> (Vec<range::Range>, HashMap<String, Vec<range::Range>>) {
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