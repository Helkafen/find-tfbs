extern crate bgzip;
extern crate clap;
extern crate which;
extern crate subprocess;
extern crate multiqueue;

use std::rc::Rc;
use rust_htslib::bcf::*;

use bgzip::write::BGzWriter;
use std::fs;
use std::fs::File;
use std::io::{BufReader,BufWriter};
use std::io::prelude::*;

use std::collections::HashMap;
use std::collections::HashSet;

use std::path::Path;
use std::time::SystemTime;

use std::sync::mpsc::{Sender, Receiver};
use std::sync::mpsc;
use std::thread;
use std::sync::{Arc,Mutex};

use clap::{Arg, App};

use subprocess::Exec;

mod types;
mod range;
mod pattern;
mod util;
mod haplotype;
mod bed;

use types::*;
use range::Range;
use pattern::{parse_pwm_files,matches};
use util::to_nucleotides_pos;
use haplotype::load_haplotypes ;
use bed::load_peak_files;


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

fn select_inner_peaks(peak: Range, peak_map: &HashMap<String, Vec<Range>>) -> HashMap<&String, Vec<&Range>> {
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

fn all_haplotype_ids(sample_count: usize) -> HashSet<HaplotypeId> {
    let mut x = HashSet::new();
    for i in (0..sample_count).into_iter() {
        x.insert(HaplotypeId {sample_id: i, side: HaplotypeSide::Left});
        x.insert(HaplotypeId {sample_id: i, side: HaplotypeSide::Right});
    }
    x
}

fn print_haplotype(haplotype: &Vec<NucleotidePos>) {
    for c in haplotype {
        print!("{}", c.nuc);
    }
    println!("");
}

fn find_all_matches(chromosome: &str, peak: &Range, reader: &mut IndexedReader, ref_haplotype: &Vec<NucleotidePos>, pwm_list: &Vec<Pattern>, mut haplotypes_with_reference_genome: HashSet<HaplotypeId>, sample_positions_in_bcf:&Vec<usize>, verbose: bool)
                    -> (Vec<Match>, u32, u32) {
    let mut match_list = Vec::new();
    let mut number_of_haplotypes = 0;
    // Find matches for people who have at least one variant inside the peak
    let (number_of_variants, mut xs) = load_haplotypes(chromosome, &peak, reader, ref_haplotype, sample_positions_in_bcf);
    if verbose {
        print!("Reference haplotype: "); print_haplotype(ref_haplotype);
    }
    for (haplotype, haplotype_ids) in xs.drain() {
        for h in haplotype_ids.iter() {
            haplotypes_with_reference_genome.remove(&h);
        }
        if verbose {
            print!("  Patched haplotype: "); print_haplotype(&haplotype);
        }
        //println!("{} PWMs", pwm_list.len());
        for pwm in pwm_list {
            match_list.extend(matches(pwm, &haplotype, haplotype_ids.clone(),verbose));
        }
        number_of_haplotypes = number_of_haplotypes + 1;
    };

    // Find matches for people who have the reference genome in this peak
    if !haplotypes_with_reference_genome.is_empty() {
        number_of_haplotypes = number_of_haplotypes + 1;
        let x: Vec<HaplotypeId> = haplotypes_with_reference_genome.into_iter().collect();
        let hap_ids = Rc::new(x);
        for pwm in pwm_list {
            match_list.extend(matches(pwm, &ref_haplotype, hap_ids.clone(),verbose));
        }
    }
    (match_list, number_of_haplotypes, number_of_variants)
}

fn read_peak_in_reference_genome(chromosome: String, peak: &Range, reference_genome: &mut bio::io::fasta::IndexedReader<std::fs::File>)-> Vec<NucleotidePos> {
    reference_genome.fetch(&chromosome, peak.start, peak.end + 1).expect("Error while seeking in reference genome file");
    let mut text = Vec::new();
    reference_genome.read(&mut text).expect(&format!("Error while reading in reference genome file {}:{}-{}", chromosome, peak.start, peak.end));
    to_nucleotides_pos(&text, &peak)
}

fn main() {
    let app = App::new("find-tfbs");
    let opt_matches: clap::ArgMatches<'static> =
                        app.version("1.0")
                        .author("Sébastian Méric de Bellefon <sebastian.meric.de.bellefon@umontreal.ca>")
                        .about("Find patterns in a VCF file")
                        .arg(Arg::with_name("chromosome")        .short("c").required(true) .takes_value(true) .value_name("CHROM")              .long("chromosome")             .help("Chromosome to scan. Ex: 'chr1'"))
                        .arg(Arg::with_name("bcf")               .short("i").required(true) .takes_value(true) .value_name("IN")                 .long("input")                  .help("BCF input file to use"))
                        .arg(Arg::with_name("output")            .short("o").required(true) .takes_value(true) .value_name("OUT")                .long("output")                 .help("Output VCF file"))
                        .arg(Arg::with_name("ref")               .short("r").required(true) .takes_value(true) .value_name("REF")                .long("reference")              .help("Reference genome. Ex: hg38.fa"))
                        .arg(Arg::with_name("bed_files")         .short("b").required(true) .takes_value(true) .value_name("BED")                .long("bed")                    .help("Bed files containing the regions to scan"))
                        .arg(Arg::with_name("pwm_names")         .short("n").required(true) .takes_value(true) .value_name("PWM_NAMES")          .long("pwm_names")              .help("List of PWM names to scan. Ex: CTCF_HUMAN.H11MO.0.A,IRF1_HUMAN.H11MO.0.A"))
                        .arg(Arg::with_name("pwm_file")          .short("p").required(true) .takes_value(true) .value_name("PWM")                .long("pwm_file")               .help("PWM file. Ex: HOCOMOCOv11_full_pwms_HUMAN_mono.txt"))
                        .arg(Arg::with_name("pwm_threshold_dir")            .required(true) .takes_value(true) .value_name("THRESHOLD_DIRECTORY").long("pwm_threshold_directory").help("PWM thresholds directory. Extracted from HOCOMOCO's thresholds_HUMAN_mono.tar.gz"))
                        .arg(Arg::with_name("pwm_threshold")     .short("t").required(true) .takes_value(true) .value_name("THRESHOLD")          .long("pwm_threshold")          .help("PWM threshold value. E.g 0.001"))
                        .arg(Arg::with_name("forward_only")      .short("f").required(false).takes_value(false)                                  .long("forward_only")           .help("Only examine the forward strand"))
                        .arg(Arg::with_name("threads")           .short("n").required(false).takes_value(true) .value_name("THREADS")            .long("threads")                .help("Size of the thread pool, in addition to the writer thread"))
                        .arg(Arg::with_name("min_maf")           .short("m").required(false).takes_value(true) .value_name("MIN_NAF")            .long("min_maf")                .help("Minimal number of occurences of the non-majority configurations"))
                        .arg(Arg::with_name("after_position")    .short("t").required(false).takes_value(true) .value_name("AFTER_POSITION")     .long("after_position")         .help("Only consider peaks that start after this position"))
                        .arg(Arg::with_name("samples")           .short("s").required(false).takes_value(true) .value_name("SAMPLES")            .long("samples")                .help("Samples file"))
                        .arg(Arg::with_name("tabix")             .short("z").required(false).takes_value(false)                                  .long("tabix")                  .help("Compress VCF with bgzip and tabix it"))
                        .arg(Arg::with_name("verbose")           .short("v").required(false).takes_value(false)                                  .long("verbose")                .help("Verbose log"))
                        .get_matches();

    let chromosome               = opt_matches.value_of("chromosome").unwrap().to_string();
    let bcf                      = opt_matches.value_of("bcf").unwrap().to_string();
    let bed_files: Vec<&str>     = opt_matches.value_of("bed_files").unwrap().split(',').collect();
    let reference_genome_file    = opt_matches.value_of("ref").unwrap().to_string();
    let pwm_file                 = opt_matches.value_of("pwm_file").unwrap();
    let pwm_threshold_directory  = opt_matches.value_of("pwm_threshold_dir").unwrap();
    let pwm_threshold: f32       = match opt_matches.value_of("pwm_threshold") {
        Some(s) => s.to_string().parse().expect("Cannot parse MAF"),
        None => panic!("No PWM threshold value"),
    };
    let wanted_pwms: Vec<String> = opt_matches.value_of("pwm_names").unwrap().split(',').map(|s| s.to_string()).collect();
    let output_file              = opt_matches.value_of("output").unwrap().to_string();
    let forward_only: bool       = opt_matches.is_present("forward_only");
    let run_tabix: bool          = opt_matches.is_present("tabix");
    let min_maf: u32             = match opt_matches.value_of("min_maf") {
        Some(s) => s.to_string().parse().expect("Cannot parse MAF"),
        None => 0,
    };

    let threads = match opt_matches.value_of("threads") {
        Some(s) => {
            let n = s.to_string().parse().expect("Cannot parse thread number");
            if n<1 { panic!("Wrong number of threads"); } else { n }
        },
        None => {1}
    };

    let after_position: u64 =
        match opt_matches.value_of("after_position") {
            Some(s) => s.to_string().parse().expect("Cannot parse after_position"),
            None => 0,
        };

    if run_tabix{
        which::which("bgzip").ok().expect("bgzip cannot in found in PATH");
        which::which("tabix").ok().expect("tabix cannot in found in PATH");
    }

    let wanted_samples = opt_matches.value_of("samples");

    let verbose: bool       = opt_matches.is_present("verbose");

    run(chromosome, bcf, bed_files, reference_genome_file, wanted_samples, pwm_file, pwm_threshold_directory, pwm_threshold, wanted_pwms, output_file, forward_only, run_tabix, min_maf, threads, after_position, verbose);

    println!("End of program.");
}

fn run(chromosome: String, bcf: String, bed_files: Vec<&str>, reference_genome_file: String, wanted_samples: Option<&str>, pwm_file: &str, pwm_threshold_directory: &str, pwm_threshold: f32,
        wanted_pwms: Vec<String>, output_file: String, forward_only: bool, run_tabix: bool, min_maf: u32, threads: u32, after_position: u64, verbose: bool) {

    let pwm_list: Vec<Pattern> = parse_pwm_files(pwm_file, pwm_threshold_directory, pwm_threshold, wanted_pwms, !forward_only);
    assert!(pwm_list.len() > 0);
    let mut pwm_name_dict: HashMap<u16, String> = HashMap::new();
    for pwm in &pwm_list {
        match pwm {
            Pattern::PWM{weights, name, pattern_id, min_score, direction} => {
                println!("PWM {} {} {} {}", name, min_score, direction, weights.len());
                pwm_name_dict.insert(*pattern_id, name.clone());
            }
            Pattern::OtherPattern{name, pattern_id} => {
                pwm_name_dict.insert(*pattern_id, name.clone());
            }
        }
    }

    let (merged_peaks, peak_map) = load_peak_files(&bed_files, &chromosome, after_position);
    let number_of_peaks = merged_peaks.len();

    let bcf_samples: Vec<String> = get_sample_names(&IndexedReader::from_path(&bcf).expect("Error while opening the bcf file"));
    let bcf_samples_len = bcf_samples.len();

    let _writer_thread = {
        // Communication channels
        let (tx, rx): (Sender<String>, Receiver<String>) = mpsc::channel();
        let (peak_tx, peak_rx) = multiqueue::broadcast_queue(merged_peaks.len() as u64);

        // The writer thread receive lines of text and writes them to the output file
        let _writer_thread = thread::spawn(move || {
            let temp = format!("{}.part", output_file.clone());
            {
                let mut writer = BGzWriter::new(BufWriter::with_capacity(4096*1000,fs::File::create(temp.clone()).expect("Could not create output file")));
                loop {
                    match rx.recv() {
                        Ok(s) => { writer.write(s.as_bytes()).expect("Could not write bytes to output file"); }
                        Err(_) => { let _ = writer.flush(); break; }
                    }
                };

                let _ = writer.flush();
            }
            if run_tabix{
                let x = Exec::shell(&format!("zcat {} | bgzip > {}; tabix -f -p vcf {}; rm {}", temp.clone(), output_file.clone(), output_file.clone(), temp.clone())).join().unwrap();
                if x.success() {
                    println!("Tabixed file {}", output_file.clone());
                }
                else {
                    println!("Failed to tabix file {}", output_file.clone());
                }
            }
            else {
                fs::rename(temp.clone(), output_file.clone()).expect(&format!("Could not rename {} into {}", temp.clone(), output_file.clone()));
            }
            println!("End of writer thread");
        });

        // Find the individuals and their column number in the BCF file
        let (samples, sample_positions_in_bcf): (Vec<String>, Vec<usize>) = {
            match wanted_samples {
            None => (bcf_samples, (0..bcf_samples_len).into_iter().collect()),
            Some(f) => {
                    let input = File::open(f).expect(&format!("Could not open sample file {}", f));
                    let buffered = BufReader::new(input);
                    let mut samples = Vec::new();
                    for line in buffered.lines() { if let Ok(l) = line { if l.len() > 1 { samples.push(l); } } }
                    let mut sample_indices = Vec::new();
                    let samples_set: HashSet<String> = samples.clone().into_iter().collect();
                    for i in 0..bcf_samples_len {
                        if samples_set.contains(&bcf_samples[i]) {
                            sample_indices.push(i);
                        }
                    }
                    (samples, sample_indices)
            },
            }
        };
        let sample_count = samples.len();
        let null_count: Vec<u32> = repeat(sample_count, 0);
        println!("Reading {} samples out of {}", sample_count, bcf_samples_len);
        let all_haplotypes_with_reference_genome: HashSet<HaplotypeId> = all_haplotype_ids(sample_count);

        // Write header in output file
        tx.send("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT".to_string()).expect("Could not create output file");
        for sample in samples {
            tx.send("\t".to_string()).expect("Could not create output file");
            tx.send(sample.clone()).expect("Could not create output file");
        }tx.send("\n".to_string()).expect("Could not create output file");

        let start_time = SystemTime::now();

        // A fake and unique position in the chromosome given for each line in the resulting vcf
        let fake_position: Arc<Mutex<u32>> = Arc::new(Mutex::new(1));
        let number_of_peaks_processed: Arc<Mutex<u32>> = Arc::new(Mutex::new(0));

        let mut handlers = vec![];
        for _ in 0..threads {
            let chromosome                = Arc::new(chromosome.clone());
            let peak_map                  = Arc::new(peak_map.clone());
            let pwm_list                  = Arc::new(pwm_list.clone());
            let number_of_peaks_processed = number_of_peaks_processed.clone();
            let sample_positions_in_bcf   = Arc::new(sample_positions_in_bcf.clone());
            let null_count                = Arc::new(null_count.clone());
            let pwm_name_dict             = Arc::new(pwm_name_dict.clone());
            let fake_position             = fake_position.clone();
            let tx                        = tx.clone();
            let peak_rx                   = peak_rx.clone();
            let all_haplotypes_with_reference_genome = Arc::new(all_haplotypes_with_reference_genome.clone());
            let mut reader = rust_htslib::bcf::IndexedReader::from_path(&bcf).expect("Error while opening the bcf file");
            let mut reference_genome = bio::io::fasta::IndexedReader::from_file(&Path::new(&reference_genome_file)).expect(&format!("Error while opening the reference genome '{}'", reference_genome_file));

            let handle = thread::spawn(move || {
                println!("Spawned one worker");
                for peaks_chunk in peak_rx {
                    for peak in peaks_chunk {
                        process_peak(&chromosome,
                            &mut reader,
                            &mut reference_genome,
                            peak,
                            &peak_map,
                            &pwm_list,
                            &sample_positions_in_bcf,
                            &null_count,
                            &pwm_name_dict,
                            min_maf,
                            &all_haplotypes_with_reference_genome,
                            start_time,
                            &tx,
                            fake_position.clone(),
                            number_of_peaks,
                            number_of_peaks_processed.clone(),
                            verbose);
                    }
                }
            });
            handlers.push(handle);
        }

        // Performance optimization: we process peaks by chunks of 50.
        // Reason: BCF files are compressed by blocks, and each block contains the genotypes for a range of coordinates
        //         Processing large contiguous chunks is faster, because we're not randomly reading lines from different compressed blocks
        let peaks_chunks: Vec<Vec<Range>> = merged_peaks.chunks(50).map(|x| x.to_vec()).collect();
        for peaks_chunk in peaks_chunks {
            let _ = peak_tx.try_send(peaks_chunk).unwrap();
        }
        drop(peak_tx); // Let the worker threads know that this is the end of the input stream

        for thread in handlers {
            let _ = thread.join();
        }
        _writer_thread
    }; // tx dies here, which allows the writer thread to terminate

    let _ = _writer_thread.join();

    println!("Writer thread joined. End program");
}

fn process_peak(chromosome: &str, reader: &mut IndexedReader, reference_genome: &mut bio::io::fasta::IndexedReader<std::fs::File>, merged_peak: Range, peak_map: &HashMap<String, Vec<range::Range>>,
                pwm_list: &Vec<Pattern>, sample_positions_in_bcf: &Vec<usize>,
                null_count: &Vec<u32>, pwm_name_dict: &HashMap<u16, String>, min_maf: u32,
                all_haplotypes_with_reference_genome: &HashSet<HaplotypeId>, start_time: SystemTime, txx: &Sender<String>,
                fake_position: Arc<Mutex<u32>>, number_of_peaks: usize, number_of_peaks_processed: Arc<Mutex<u32>>, verbose: bool) -> () {
    let peak_start_time = SystemTime::now();

    let chr = chromosome.replace("chr", "");

    let largest_pwm_size: u32 = pwm_list.iter().map(|x| pattern_length(x.clone())).max().unwrap();

    // The TFBS may overlap the right and left border of the merged peak, so we build larger haplotypes
    let extended_merged_peak = range::Range::new(merged_peak.start-largest_pwm_size as u64 + 1,merged_peak.end+largest_pwm_size as u64 - 1);

    let ref_haplotype = read_peak_in_reference_genome(chromosome.to_string(), &extended_merged_peak, reference_genome);

    let inner_peaks: HashMap<&String, Vec<&Range>> = select_inner_peaks(merged_peak, &peak_map);

    let (match_list, number_of_haplotypes, number_of_variants) = find_all_matches(&chromosome, &extended_merged_peak, reader, &ref_haplotype, &*pwm_list, (*all_haplotypes_with_reference_genome).clone(), &*sample_positions_in_bcf, verbose);

    for ((source, inner_peak, pattern_id),(v1,v2)) in count_matches_by_sample(&match_list, &inner_peaks, &*null_count).drain() {
        let pwm_name = pwm_name_dict.get(&pattern_id).expect("Logic error: No pattern name for a pattern_id");
        let id_str = format!("{},{},{}-{}",source, pwm_name, inner_peak.start, inner_peak.end);

        //println!("count_matches_by_sample returned something {}", id_str);
        if let Some((distinct_counts, maf, freq0, freq1, freq2, genotypes)) = counts_as_genotypes(v1, v2, verbose) {
            if maf >= min_maf {
                let distinct_counts_str: Vec<String> = distinct_counts.iter().map(|c| c.to_string()).collect();
                let info_str = format!("COUNTS={};freqs={}/{}/{}", distinct_counts_str.join(","), freq0, freq1, freq2);
                txx.send(format!("{}\t{}\t{}\t.\t.\t.\t.\t{}\tGT:DS{}\n", chr, fake_position.lock().unwrap(), id_str, info_str, genotypes).to_string()).expect("Could not write result");
                *fake_position.lock().unwrap() += 1;
            }
            else { if verbose { println!("Frequency insufficient"); } }
        }
    }

    let number_of_matches: u64 = match_list.iter().map(|m| m.haplotype_ids.len() as u64).sum();
    let peak_time_elapsed = peak_start_time.elapsed().unwrap().as_millis();
    let global_time_elapsed = start_time.elapsed().unwrap().as_millis();
    *number_of_peaks_processed.lock().unwrap() += 1;
    println!("Peak {}/{}\t{} ms ({} total)\t{}\t{}\t{} haplotypes\t{} variants\t{} matches", number_of_peaks_processed.lock().unwrap(), number_of_peaks, peak_time_elapsed, global_time_elapsed, merged_peak.start, merged_peak.end, number_of_haplotypes, number_of_variants, number_of_matches);
}


fn counts_as_genotypes(v1: Vec<u32>, v2: Vec<u32>, verbose: bool) -> Option<(Vec<u32>, u32, u32, u32, u32, String)> {
    // Sum of the number of TFBS on both alleles
    let v = {
        assert!(v1.len() == v2.len());
        let mut v = v1.clone();
        for (i,x) in v2.iter().enumerate() {
            v[i] += x;
        }
        v
    };

    let min = v.iter().min();
    let max = v.iter().max();

    match (min, max) {
        (Some(&lowest), Some(&highest)) => {
            if verbose { println!("Min and max count: {} {}", lowest, highest); }
            if lowest == highest {
                None // No variation in the number of TFBS, so let's forget about this region
            }
            else {
                let mut res = String::with_capacity(v.len()*4);
                let intermediate_1_1000 = (lowest * 1000 * 3 + highest * 1000    ) / 4;
                let intermediate_3_1000 = (lowest * 1000     + highest * 1000 * 3) / 4;
                let mut all_values = vec![lowest, highest];
                let mut zero_count: u32 = 0;
                let mut one_count: u32 = 0;
                let mut two_count: u32 = 0;
                let lowest_f32 = lowest as f32;
                let spread_f32 = highest as f32 - lowest_f32;
                for x in v {
                    if x == lowest       { res.push_str("\t0|0:0.0"); zero_count += 1; }
                    else if x == highest { res.push_str("\t1|1:2.0"); two_count += 1;  }
                    else {
                        if !all_values.contains(&x) { all_values.push(x);}
                        let x_1000 = x *1000;
                        if x_1000 < intermediate_1_1000 { res.push_str("\t0|0"); zero_count += 1; }
                        else if x_1000 < intermediate_3_1000 { res.push_str("\t0|1"); one_count += 1; }
                        else { res.push_str("\t1|1"); two_count += 1; }
                        let pseudo_dosage = ((x as f32 - lowest_f32)*2.0)/ spread_f32; // Simulating a dosage, between 0 and 2. EMMAX accepts this field
                        res.push_str(&format!(":{:.4}", pseudo_dosage));
                    }
                }
                let maf =
                    if zero_count >= one_count && zero_count >= two_count {
                        one_count + two_count
                    }
                    else if two_count >= zero_count && two_count >= one_count {
                        zero_count + one_count
                    }
                    else { zero_count + two_count };
                all_values.sort();
                Some((all_values, maf, zero_count, one_count, two_count, res))
            }
        },
        (None, _) => None,
        (_, None) => None,
    }

}

fn count_matches_by_sample<'a>(match_list: &Vec<Match>, inner_peaks: &'a HashMap<&String, Vec<&Range>>, null_count: &Vec<u32>) -> HashMap<(&'a String, &'a Range, u16), (Vec<u32>,Vec<u32>)> {
    let mut pppp: HashMap<(&String, &Range, u16), (Vec<u32>, Vec<u32>)> = HashMap::new();
    for m in match_list {
        for (&s,inner) in inner_peaks.iter().map(|(s,x)| (s, x.iter().filter(|y| y.overlaps(&m.range)))) {
            for &inner_peak in inner {
                let key = (s, inner_peak, m.pattern_id);
                match pppp.get_mut(&key) {
                    Some((l,r)) => {
                        for haplotype_id in m.haplotype_ids.iter() {
                            if haplotype_id.side == HaplotypeSide::Left {
                                l[haplotype_id.sample_id as usize] = l[haplotype_id.sample_id as usize] + 1;
                            }
                            else {
                                r[haplotype_id.sample_id as usize] = r[haplotype_id.sample_id as usize] + 1;
                            }
                        }
                    }
                    None => {
                        let (mut l, mut r) = (null_count.clone(), null_count.clone());
                        for haplotype_id in m.haplotype_ids.iter() {
                            if haplotype_id.side == HaplotypeSide::Left {
                                l[haplotype_id.sample_id as usize] = l[haplotype_id.sample_id as usize] + 1;
                            }
                            else {
                                r[haplotype_id.sample_id as usize] = r[haplotype_id.sample_id as usize] + 1;
                            }
                        }
                        pppp.insert(key, (l,r));
                    }
                }
            }
        }
    }
    pppp
}


#[cfg(test)]
mod tests {
    use super::*;

    fn wd() -> String {
        let current_dir = std::env::current_dir().unwrap();
        let current_dir_str = current_dir.as_path().to_str().unwrap();
        current_dir_str.to_string()
    }


    #[test]
    fn test_integration_no_polymorphism() {
        run("chr1".to_string(), "test_data/genotypes.bcf".to_string(), vec!["test_data/regions1.bed","test_data/regions2.bed"], "test_data/reference_genome.fa".to_string(), 
             Some("test_data/samples"),
             "test_data/pwm_definitions.txt",
             "test_data", 0.0001, vec!["ACGT".to_string()], "test_data/output1.vcf.gz".to_string(), false, false, 0, 1, 0, true);
        let output = std::fs::read(format!("{}/test_data/output1.vcf.gz", wd())).expect("unable to read file");
        let expected = std::fs::read(format!("{}/test_data/expected_output_1.vcf.gz", wd())).expect("unable to read file");
        assert_eq!(output, expected);
    }

    #[test]
    fn test_integration_one_polymorphism() {
        run("chr1".to_string(), "test_data/genotypes2.bcf".to_string(), vec!["test_data/regions1.bed","test_data/regions2.bed"], "test_data/reference_genome.fa".to_string(), 
             Some("test_data/samples"),
             "test_data/pwm_definitions.txt",
             "test_data", 0.0001, vec!["ACGT".to_string()], "test_data/output2.vcf.gz".to_string(), false, false, 0, 1, 0, true);
        let output = std::fs::read(format!("{}/test_data/output2.vcf.gz", wd())).expect("unable to read file");
        let expected = std::fs::read(format!("{}/test_data/expected_output_2.vcf.gz", wd())).expect("unable to read file");
        assert_eq!(output, expected);
    }

    #[test]
    fn test_count_matches() {
        let sample_count = 2;

        let match1 = Match { range: range::Range::new(10,11), pattern_id: 0, haplotype_ids: Rc::new(vec![HaplotypeId { sample_id: 0, side: HaplotypeSide::Left } ]) };

        let match_list1 = vec![match1.clone()];
        let match_list2 = vec![Match{range: range::Range::new(20,21), ..match1.clone()}];
        let match_list3 = vec![Match{range: range::Range::new(4, 5) , ..match1.clone()}];
        let match_list4 = vec![Match{range: range::Range::new(3, 4) , ..match1.clone()}];
        let match_list5 = vec![Match{range: range::Range::new(21, 22) , ..match1.clone()}];

        let match_list6 = vec![
            Match { range: range::Range::new(4,5), pattern_id: 9, haplotype_ids: Rc::new(vec![HaplotypeId { sample_id: 1, side: HaplotypeSide::Right } ]) } 
        ];

        let match_list7 = vec![
            Match { range: range::Range::new(17,18), pattern_id: 11, haplotype_ids: Rc::new(vec![HaplotypeId { sample_id: 1, side: HaplotypeSide::Right } ]) } 
        ];

        let mep = "MEP".to_string();
        let erythro = "Erythro".to_string();
        let range1 = range::Range::new(5,20);
        let range2 = range::Range::new(15,25);

        let inner_peaks = {
            let mut x = HashMap::new();
            x.insert(&mep, vec![&range1]);
            x
        };

        let inner_peaks2 = {
            let mut x = HashMap::new();
            x.insert(&mep, vec![&range1]);
            x.insert(&erythro, vec![&range2]);
            x
        };

        let null_count: Vec<u32> = repeat(sample_count, 0);

        // TFBS overlaps open chromatin region
        assert_eq!(count_matches_by_sample(&match_list1, &inner_peaks, &null_count), 
            {
                let mut x = HashMap::new();
                x.insert((&mep, &range1, 0), (vec![1,0], vec![0,0]));
                x
            });
        assert_eq!(count_matches_by_sample(&match_list1, &inner_peaks, &null_count), count_matches_by_sample(&match_list2, &inner_peaks, &null_count));
        assert_eq!(count_matches_by_sample(&match_list1, &inner_peaks, &null_count), count_matches_by_sample(&match_list3, &inner_peaks, &null_count));

        // TFBS doesn't overlap open chromatin region
        assert_eq!(count_matches_by_sample(&match_list4, &inner_peaks, &null_count), count_matches_by_sample(&match_list5, &inner_peaks, &null_count));
        assert_eq!(count_matches_by_sample(&match_list4, &inner_peaks, &null_count), HashMap::new());

        // Use a different pattern_id and move the TFBS to another person/chromosome
        assert_eq!(count_matches_by_sample(&match_list6, &inner_peaks, &null_count), 
            {
                let mut x = HashMap::new();
                x.insert((&mep, &range1, 9), (vec![0,0], vec![0,1]));
                x
            });

        assert_eq!(count_matches_by_sample(&match_list1, &inner_peaks2, &null_count), 
        {
            let mut x = HashMap::new();
            x.insert((&mep, &range1, 0), (vec![1,0], vec![0,0]));
            x
        });

        assert_eq!(count_matches_by_sample(&match_list2, &inner_peaks2, &null_count),
        {
            let mut x = HashMap::new();
            x.insert((&mep, &range1, 0), (vec![1,0], vec![0,0]));
            x.insert((&erythro, &range2, 0), (vec![1,0], vec![0,0]));
            x
        });
        assert_eq!(count_matches_by_sample(&match_list1, &inner_peaks2, &null_count), count_matches_by_sample(&match_list3, &inner_peaks2, &null_count));
        assert_eq!(count_matches_by_sample(&match_list4, &inner_peaks2, &null_count), HashMap::new());
        assert_eq!(count_matches_by_sample(&match_list5, &inner_peaks2, &null_count),
        {
            let mut x = HashMap::new();
            x.insert((&erythro, &range2, 0), (vec![1,0], vec![0,0]));
            x
        });
        assert_eq!(count_matches_by_sample(&match_list4, &inner_peaks2, &null_count), HashMap::new());

        // Use a different pattern_id and move the TFBS to another person/chromosome
        assert_eq!(count_matches_by_sample(&match_list6, &inner_peaks2, &null_count), 
        {
            let mut x = HashMap::new();
            x.insert((&mep, &range1, 9), (vec![0,0], vec![0,1]));
            x
        });

        assert_eq!(count_matches_by_sample(&match_list7, &inner_peaks2, &null_count), 
        {
            let mut x = HashMap::new();
            x.insert((&mep, &range1, 11), (vec![0,0], vec![0,1]));
            x.insert((&erythro, &range2, 11), (vec![0,0], vec![0,1]));
            x
        });
    }

}
