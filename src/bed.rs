use super::range;

use bio::io::bed;
use std::collections::HashMap;
use std::path::Path;



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

pub fn load_peak_files(bed_files: &Vec<&str>, chromosome: &str, after_position: u64) -> (Vec<range::Range>, HashMap<String, Vec<range::Range>>) {
    let mut peak_map: HashMap<String, Vec<range::Range>> = HashMap::new();
    for bed_file in bed_files {
        if Path::new(bed_file).exists() {
            let peaks = load_bed(&bed_file.to_string(), chromosome);
            println!("Loaded {}:\t {} peaks covering {} bp", bed_file, peaks.len(), sum_peak_sizes(&peaks));
            peak_map.insert(bed_file.to_string(), peaks.into_iter().filter(|p| p.start >= after_position).collect());
        }
        else {
            panic!("Bed file {} does not exist", bed_file);
        }
    }
    let vals: Vec<Vec<range::Range>> = peak_map.values().map(|x| x.clone()).collect();
    let rs: range::RangeStack = vals.concat().iter().collect();
    println!("Merged all peak files: {} peaks covering {} bp", rs.ranges.len(), sum_peak_sizes(&rs.ranges));
    (rs.ranges, simplify_peak_map(peak_map))
}

fn simplify(s: String) -> String {
    let p = Path::new(&s).file_name().unwrap();
    String::from(p.to_str().unwrap())
}

fn simplify_peak_map(mut m: HashMap<String, Vec<range::Range>>) -> HashMap<String, Vec<range::Range>> {
    let mut res = HashMap::new();
    for (k,v) in m.drain() {
        res.insert(simplify(k), v);
    }
    res
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merge_bed() {
        let chromosome = "chr1";
        let beds = vec!["test_data/regions1.bed", "test_data/regions2.bed"];
        let (merged_ranges, range_map) = load_peak_files(&beds, &chromosome, 0);

        assert_eq!(merged_ranges[0], range::Range::new(100, 115));
        assert_eq!(merged_ranges[1], range::Range::new(118, 130));
        assert_eq!(merged_ranges[2], range::Range::new(150, 160));
        assert_eq!(merged_ranges[3], range::Range::new(161, 165));
        assert_eq!(merged_ranges[4], range::Range::new(180, 210));

        let expected = {
            let mut ex = HashMap::new();
            let region_1 = vec![range::Range::new(100,110),
                                range::Range::new(120,130),
                                range::Range::new(150,160),
                                range::Range::new(180,190),
                                range::Range::new(200,210)];
            let region_2 = vec![range::Range::new(110,115),
                                range::Range::new(118,125),
                                range::Range::new(161,165),
                                range::Range::new(190,200)];
            ex.insert("regions1.bed".to_string(), region_1);
            ex.insert("regions2.bed".to_string(), region_2);
            ex
        };
        assert_eq!(range_map, expected);
        assert_eq!(sum_peak_sizes(&merged_ranges), 71);
    }
}