extern crate bgzip;
extern crate clap;
extern crate rayon;
extern crate which;
extern crate subprocess;

use std::rc::Rc;
use rust_htslib::bcf::*;

use bgzip::write::BGzWriter;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

use std::collections::HashMap;
use std::collections::HashSet;

use std::path::Path;
use std::time::SystemTime;

use std::sync::mpsc::{Sender, Receiver};
use std::sync::mpsc;
use std::thread;
use std::sync::{Arc,Mutex};
use rayon::prelude::*;

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
use haplotype::load_haplotypes;
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

fn find_all_matches(chromosome: &str, peak: &Range, reader: &mut IndexedReader, ref_haplotype: &Vec<NucleotidePos>, pwm_list: &Vec<PWM>, mut haplotypes_with_reference_genome: HashSet<HaplotypeId>, sample_positions_in_bcf:&Vec<usize>)
                    -> (Vec<Match>, u32, u32) {
    let mut match_list = Vec::new();
    let mut number_of_haplotypes = 0;
    // Find matches for people who have at least one variant inside the peak
    let (number_of_variants, mut xs) = load_haplotypes(chromosome, &peak, reader, ref_haplotype, sample_positions_in_bcf);
    for (haplotype, haplotype_ids) in xs.drain() {
        for h in haplotype_ids.iter() {
            haplotypes_with_reference_genome.remove(&h);
        }
        for pwm in pwm_list {
            match_list.extend(matches(pwm, &haplotype, haplotype_ids.clone()));
        }
        number_of_haplotypes = number_of_haplotypes + 1;
    };

    // Find matches for people who have the reference genome in this peak
    if !haplotypes_with_reference_genome.is_empty() {
        number_of_haplotypes = number_of_haplotypes + 1;
        let x: Vec<HaplotypeId> = haplotypes_with_reference_genome.into_iter().collect();
        let hap_ids = Rc::new(x);
        for pwm in pwm_list {
            match_list.extend(matches(pwm, &ref_haplotype, hap_ids.clone()));
        }
    }
    (match_list, number_of_haplotypes, number_of_variants)
}

fn read_peak_in_reference_genome(chromosome: &str, peak: &Range, reference_genome: &mut bio::io::fasta::IndexedReader<std::fs::File>)-> Vec<NucleotidePos> {
    reference_genome.fetch(chromosome, peak.start, peak.end - 1).expect("Error while seeking in reference genome file");
    let mut text = Vec::new();
    reference_genome.read(&mut text).expect(&format!("Error while reading in reference genome file {}:{}-{}", chromosome, peak.start, peak.end));
    to_nucleotides_pos(&text, &peak)
}

// cargo run --release -- -c chr1 -m 5 -i /home/seb/masters/topmed/source/TOPMed_dbGaP_20180710/dbGaP-12336/68779/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.8/phased/freeze.8.chr1.pass_only.phased.bcf -b bed/Bcell-13.bed,bed/CD4-9.bed,bed/CD8-10.bed,bed/CLP-14.bed,bed/CMP-4.bed,bed/Erythro-15.bed,bed/GMP-5.bed,bed/HSC-1.bed,bed/LMPP-3.bed,bed/MCP.bed,bed/mDC.bed,bed/MEGA1.bed,bed/MEGA2.bed,bed/MEP-6.bed,bed/Mono-7.bed,bed/MPP-2.bed,bed/Nkcell-11.bed,bed/pDC.bed -o test2.gz -p /home/seb/masters/regu/dnamotifs/HOCOMOCOv11_full_pwms_HUMAN_mono.txt -n JUNB_HUMAN.H11MO.0.A,FOSL1_HUMAN.H11MO.0.A,FOSL2_HUMAN.H11MO.0.A,JDP2_HUMAN.H11MO.0.D,GATA1_HUMAN.H11MO.0.A,GATA2_HUMAN.H11MO.0.A,GATA3_HUMAN.H11MO.0.A,GATA4_HUMAN.H11MO.0.A,GATA5_HUMAN.H11MO.0.D,GATA6_HUMAN.H11MO.0.A,JUN_HUMAN.H11MO.0.A,JUND_HUMAN.H11MO.0.A,BATF_HUMAN.H11MO.0.A,ATF3_HUMAN.H11MO.0.A,BACH1_HUMAN.H11MO.0.A,BACH2_HUMAN.H11MO.0.A,NFE2_HUMAN.H11MO.0.A,CEBPA_HUMAN.H11MO.0.A,CEBPB_HUMAN.H11MO.0.A,CEBPD_HUMAN.H11MO.0.C,CEBPE_HUMAN.H11MO.0.A,CEBPG_HUMAN.H11MO.0.B,SPIB_HUMAN.H11MO.0.A,IRF8_HUMAN.H11MO.0.B,SPI1_HUMAN.H11MO.0.A,MESP1_HUMAN.H11MO.0.D,ID4_HUMAN.H11MO.0.D,HTF4_HUMAN.H11MO.0.A,ITF2_HUMAN.H11MO.0.C,STAT1_HUMAN.H11MO.0.A,STAT2_HUMAN.H11MO.0.A,SPIC_HUMAN.H11MO.0.D,CTCF_HUMAN.H11MO.0.A,IRF1_HUMAN.H11MO.0.A,DBP_HUMAN.H11MO.0.B,MAFK_HUMAN.H11MO.1.A,ATF4_HUMAN.H11MO.0.A,ASCL1_HUMAN.H11MO.0.A,ASCL2_HUMAN.H11MO.0.D,TFE2_HUMAN.H11MO.0.A,MYOD1_HUMAN.H11MO.0.A,EVI1_HUMAN.H11MO.0.B,IRF3_HUMAN.H11MO.0.B,ZEB1_HUMAN.H11MO.0.A,IRF9_HUMAN.H11MO.0.C,HEN1_HUMAN.H11MO.0.C,LYL1_HUMAN.H11MO.0.A -t /home/seb/masters/regu/dnamotifs/hocomoco_thresholds.tab -r /home/seb/masters/hg38.fa

fn main() {
    let opt_matches = App::new("VCF_PWM")
                        .version("1.0")
                        .author("Sébastian Méric de Bellefon <arnaudpourseb@gmail.com>")
                        .about("Find patterns in a VCF file")
                        .arg(Arg::with_name("chromosome")        .short("c").required(true) .takes_value(true) .value_name("CHROM")       .long("chromosome")        .help("Chromosome to scan. Ex: 'chr1'"))
                        .arg(Arg::with_name("bcf")               .short("i").required(true) .takes_value(true) .value_name("IN")          .long("input")             .help("BCF input file to use"))
                        .arg(Arg::with_name("output")            .short("o").required(true) .takes_value(true) .value_name("OUT")         .long("output")            .help("Output VCF file"))
                        .arg(Arg::with_name("ref")               .short("r").required(true) .takes_value(true) .value_name("REF")         .long("reference")         .help("Reference genome. Ex: hg38.fa"))
                        .arg(Arg::with_name("bed_files")         .short("b").required(true) .takes_value(true) .value_name("BED")         .long("bed")               .help("Bed files containing the regions to scan"))
                        .arg(Arg::with_name("pwm_names")         .short("n").required(true) .takes_value(true) .value_name("PWM_NAMES")   .long("pwm_names")         .help("List of PWM names to scan. Ex: CTCF_HUMAN.H11MO.0.A,IRF1_HUMAN.H11MO.0.A"))
                        .arg(Arg::with_name("pwm_file")          .short("p").required(true) .takes_value(true) .value_name("PWM")         .long("pwm_file")          .help("PWM file. Ex: HOCOMOCOv11_full_pwms_HUMAN_mono.txt"))
                        .arg(Arg::with_name("pwm_threshold_file").short("t").required(true) .takes_value(true) .value_name("THRESHOLD")   .long("pwm_threshold_file").help("PWM threshold file. Ex: HOCOMOCOv11_full_standard_thresholds_HUMAN_mono.txt"))
                        .arg(Arg::with_name("forward_only")      .short("f").required(false).takes_value(false)                           .long("forward_only")      .help("Only examine the forward strand"))
                        .arg(Arg::with_name("threads")           .short("n").required(false).takes_value(true) .value_name("THREADS")     .long("threads")           .help("Size of the thread pool, in addition to the writer thread"))
                        .arg(Arg::with_name("min_maf")           .short("m").required(false).takes_value(true) .value_name("MIN_NAF")     .long("min_maf")           .help("Minimal number of occurences of the non-majority configurations"))
                        .arg(Arg::with_name("samples")           .short("s").required(false).takes_value(true) .value_name("SAMPLES")     .long("samples")           .help("Samples file"))
                        .arg(Arg::with_name("tabix")             .short("z").required(false).takes_value(false)                           .long("tabix")             .help("Compress VCF with bgzip and tabix it"))
                        .get_matches();

    let chromosome               = opt_matches.value_of("chromosome").unwrap();                     //1
    let bcf                      = opt_matches.value_of("bcf").unwrap();                            //format!("/home/seb/masters/topmed/source/TOPMed_dbGaP_20180710/dbGaP-12336/65066/topmed-dcc/exchange/phs000964_TOPMed_WGS_JHS/Combined_Study_Data/Genotypes/freeze.6a/phased/freeze.6a.{}.pass_only.phased.bcf", chromosome);
    let bed_files: Vec<&str>     = opt_matches.value_of("bed_files").unwrap().split(',').collect(); //bed/Bcell-13.bed,bed/CD4-9.bed,bed/CD8-10.bed,bed/CLP-14.bed,bed/CMP-4.bed,bed/Erythro-15.bed,bed/GMP-5.bed,bed/HSC-1.bed,bed/LMPP-3.bed,bed/MCP.bed,bed/mDC.bed,bed/MEGA1.bed,bed/MEGA2.bed,bed/MEP-6.bed,bed/Mono-7.bed,bed/MPP-2.bed,bed/Nkcell-11.bed,bed/pDC.bed
    let reference_genome_file    = opt_matches.value_of("ref").unwrap();                            //"/home/seb/masters/hg38.fa";
    let pwm_file                 = opt_matches.value_of("pwm_file").unwrap();                       //"/home/seb/masters/regu/dnamotifs/HOCOMOCOv11_full_pwms_HUMAN_mono.txt";
    let pwm_threshold_file       = opt_matches.value_of("pwm_threshold_file").unwrap();             //"/home/seb/masters/regu/dnamotifs/hocomoco_thresholds.tab";
    let wanted_pwms: Vec<String> = opt_matches.value_of("pwm_names").unwrap().split(',').map(|s| s.to_string()).collect(); //"JUNB_HUMAN.H11MO.0.A,FOSL1_HUMAN.H11MO.0.A,FOSL2_HUMAN.H11MO.0.A,JDP2_HUMAN.H11MO.0.D,GATA1_HUMAN.H11MO.0.A,GATA2_HUMAN.H11MO.0.A,GATA3_HUMAN.H11MO.0.A,GATA4_HUMAN.H11MO.0.A,GATA5_HUMAN.H11MO.0.D,GATA6_HUMAN.H11MO.0.A,JUN_HUMAN.H11MO.0.A,JUND_HUMAN.H11MO.0.A,BATF_HUMAN.H11MO.0.A,ATF3_HUMAN.H11MO.0.A,BACH1_HUMAN.H11MO.0.A,BACH2_HUMAN.H11MO.0.A,NFE2_HUMAN.H11MO.0.A,CEBPA_HUMAN.H11MO.0.A,CEBPB_HUMAN.H11MO.0.A,CEBPD_HUMAN.H11MO.0.C,CEBPE_HUMAN.H11MO.0.A,CEBPG_HUMAN.H11MO.0.B,SPIB_HUMAN.H11MO.0.A,IRF8_HUMAN.H11MO.0.B,SPI1_HUMAN.H11MO.0.A,MESP1_HUMAN.H11MO.0.D,ID4_HUMAN.H11MO.0.D,HTF4_HUMAN.H11MO.0.A,ITF2_HUMAN.H11MO.0.C,STAT1_HUMAN.H11MO.0.A,STAT2_HUMAN.H11MO.0.A,SPIC_HUMAN.H11MO.0.D,CTCF_HUMAN.H11MO.0.A,IRF1_HUMAN.H11MO.0.A,DBP_HUMAN.H11MO.0.B,MAFK_HUMAN.H11MO.1.A,ATF4_HUMAN.H11MO.0.A,ASCL1_HUMAN.H11MO.0.A,ASCL2_HUMAN.H11MO.0.D,TFE2_HUMAN.H11MO.0.A,MYOD1_HUMAN.H11MO.0.A,EVI1_HUMAN.H11MO.0.B,IRF3_HUMAN.H11MO.0.B,ZEB1_HUMAN.H11MO.0.A,IRF9_HUMAN.H11MO.0.C,HEN1_HUMAN.H11MO.0.C,LYL1_HUMAN.H11MO.0.A".split(',').into_iter().map(|a| a.to_string()).collect();
    let output_file              = opt_matches.value_of("output").unwrap().to_string();             //"test2.gz";
    let forward_only: bool       = opt_matches.is_present("forward_only");
    let run_tabix: bool          = opt_matches.is_present("tabix");
    let min_maf: u32             = match opt_matches.value_of("min_maf") {
        Some(s) => s.to_string().parse().expect("Cannot parse MAF"),
        None => 0,
    };

    if let Some(s) = opt_matches.value_of("threads") {
        let n = s.to_string().parse().expect("Cannot parse thread number");
        if n<1 { panic!("Wrong number of threads"); } else { rayon::ThreadPoolBuilder::new().num_threads(n).build_global().expect("Couldn't build the thread pool"); }
    }

    if run_tabix{
        which::which("bgzip").ok().expect("bgzip cannot in found in PATH");
        which::which("tabix").ok().expect("tabix cannot in found in PATH");
    }

    let pwm_list: Vec<PWM> = parse_pwm_files(pwm_file, pwm_threshold_file, wanted_pwms, !forward_only);
    let mut pwm_name_dict = HashMap::new();
    for pwm in &pwm_list {
        println!("PWM {} {} {}", pwm.name, pwm.min_score, pwm.direction);
        pwm_name_dict.insert(pwm.pattern_id, (pwm.name.clone(), pwm.direction.clone()));
    }

    let (merged_peaks, peak_map) = load_peak_files(&bed_files, chromosome);

    let (tx, rx): (Sender<String>, Receiver<String>) = mpsc::channel();

    let _writer_thread = thread::spawn(move || {
        let mut writer = BGzWriter::new(fs::File::create(output_file.clone() + ".part").expect("Could not create output file"));
        loop {
            match rx.recv() {
                Ok(s) => { writer.write(s.as_bytes()).expect("Could not write bytes to output file"); }
                Err(_) => { break; }
            }
        };
        if run_tabix{
            println!("Tabixing {}", output_file.clone());
            let bgzip_exit_status = Exec::shell(&format!("zcat {}.part | bgzip > {}; tabix -f -p vcf {}; rm {}.part", output_file.clone(), output_file.clone(), output_file.clone(), output_file.clone())).join().unwrap();
            println!("HERE2 {:#?}", bgzip_exit_status); // TODO: Mysteriously, we never reach this point
        }
        else {
            fs::rename(output_file.clone() + ".part", output_file.clone()).expect(&format!("Count not rename {} into {}", output_file.clone(), output_file.clone() + ".part"));
        }
    });

    let mut reader = IndexedReader::from_path(bcf).expect("Error while opening the bcf file");
    let bcf_samples = get_sample_names(&mut reader);
    let bcf_samples_length = bcf_samples.len();
    let number_of_peaks = &merged_peaks.len();

    let (samples, sample_positions_in_bcf): (Vec<String>, Vec<usize>) = {
        match opt_matches.value_of("samples") {
           None => (bcf_samples, (0..bcf_samples_length).into_iter().collect()),
           Some(f) => {
                let input = File::open(f).expect(&format!("Could not open sample file {}", f));
                let buffered = BufReader::new(input);
                let mut samples = Vec::new();
                for line in buffered.lines() { if let Ok(l) = line { if l.len() > 1 { samples.push(l); } } }
                let mut sample_indices = Vec::new();
                let samples_set: HashSet<String> = samples.clone().into_iter().collect();
                for i in 0..bcf_samples_length {
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
    println!("Reading {} samples out of {}", sample_count, bcf_samples_length);
    let all_haplotypes_with_reference_genome: HashSet<HaplotypeId> = all_haplotype_ids(sample_count);

    // Write header in output file
    tx.send("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT".to_string()).expect("Could not create output file");
    for sample in samples {
        tx.send("\t".to_string()).expect("Could not create output file");
        tx.send(sample.clone()).expect("Could not create output file");
    }tx.send("\n".to_string()).expect("Could not create output file");

    // A fake and unique position in the chromosome given for each line in the resulting vcf
    let fake_position: Arc<Mutex<u32>> = Arc::new(Mutex::new(1));
    let number_of_peaks_processed: Arc<Mutex<u32>> = Arc::new(Mutex::new(0));
    let chr = String::from(chromosome).replace("chr", "");

    let start_time = SystemTime::now();

    merged_peaks.into_par_iter().for_each_with(tx, |txx, peak| {
        let mut reader = IndexedReader::from_path(bcf).expect("Error while opening the bcf file");
        let mut reference_genome = bio::io::fasta::IndexedReader::from_file(&Path::new(reference_genome_file)).expect("Error while opening the reference genome");

        let peak_start_time = SystemTime::now();

        let ref_haplotype = read_peak_in_reference_genome(chromosome, &peak, &mut reference_genome);

        let inner_peaks: HashMap<&String, Vec<&Range>> = select_inner_peaks(peak, &peak_map);

        let (match_list, number_of_haplotypes, number_of_variants) = find_all_matches(chromosome, &peak, &mut reader, &ref_haplotype, &pwm_list, all_haplotypes_with_reference_genome.clone(), &sample_positions_in_bcf);

        for ((source, inner_peak, pattern_id),v) in count_matches_by_sample(&match_list, &inner_peaks, &null_count).drain() {
            let (distinct_counts, maf, freq0, freq1, freq2, genotypes) = counts_as_genotypes(v);
            if maf >= min_maf && distinct_counts.len() > 1 {
                let (pwm_name, pwm_direction) = pwm_name_dict.get(&pattern_id).expect("Logic error: No pattern name for a pattern_id");
                let id_str = format!("{},{},{},{}-{}",source, pwm_name, pwm_direction, inner_peak.start, inner_peak.end);
                let distinct_counts_str: Vec<String> = distinct_counts.iter().map(|c| c.to_string()).collect();
                let info_str = format!("COUNTS={};freqs={}/{}/{}", distinct_counts_str.join(","), freq0, freq1, freq2);
                txx.send(format!("{}\t{}\t{}\t.\t.\t.\t.\t{}\tGT{}\n", chr, fake_position.lock().unwrap(), id_str, info_str, genotypes).to_string()).expect("Could not write result");
                *fake_position.lock().unwrap() += 1;
            }
        }

        let number_of_matches: u64 = match_list.iter().map(|m| m.haplotype_ids.len() as u64).sum();
        let peak_time_elapsed = peak_start_time.elapsed().unwrap().as_millis();
        let global_time_elapsed = start_time.elapsed().unwrap().as_millis();
        *number_of_peaks_processed.lock().unwrap() += 1;
        println!("Peak {}/{}\t{} ms ({} total)\t{}\t{}\t{} haplotypes\t{} variants\t{} matches", number_of_peaks_processed.lock().unwrap(), number_of_peaks, peak_time_elapsed, global_time_elapsed, peak.start, peak.end, number_of_haplotypes, number_of_variants, number_of_matches);
    });
}


fn counts_as_genotypes(v: Vec<u32>) -> (Vec<u32>, u32, u32, u32, u32, String) {
    let mut res = String::with_capacity(v.len()*4);
    let min = v.iter().min();
    let max = v.iter().max();
    match (min, max) {
        (Some(&lowest), Some(&highest)) => {
            let intermediate_1_1000 = (lowest * 1000 * 3 + highest * 1000    ) / 4;
            let intermediate_3_1000 = (lowest * 1000     + highest * 1000 * 3) / 4;
            let mut all_values = vec![lowest, highest];
            let mut zero_count: u32 = 0;
            let mut one_count: u32 = 0;
            let mut two_count: u32 = 0;
            for x in v {
                if x == lowest {res.push_str("\t0|0"); zero_count += 1; }
                else if x == highest {res.push_str("\t1|1"); two_count += 1; }
                else {
                    if !all_values.contains(&x) { all_values.push(x);}
                    let x_1000 = x *1000;
                    if x_1000 < intermediate_1_1000 { res.push_str("\t0|0"); zero_count += 1; }
                    else if x_1000 < intermediate_3_1000 { res.push_str("\t0|1"); one_count += 1; }
                    else { res.push_str("\t1|1"); two_count += 1; }
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
            (all_values, maf, zero_count, one_count, two_count, res)
        },
        (None, _) => (Vec::new(), 0, 0, 0, 0, String::new()),
        (_, None) => (Vec::new(), 0, 0, 0, 0, String::new()),
    }
}

fn count_matches_by_sample<'a>(match_list: &Vec<Match>, inner_peaks: &'a HashMap<&String, Vec<&Range>>, null_count: &Vec<u32>) -> HashMap<(&'a String, &'a Range, u16), Vec<u32>> {
    let mut pppp: HashMap<(&String, &Range, u16), Vec<u32>> = HashMap::new();
    for m in match_list {
        let pos = m.pos;
        for (&s,inner) in inner_peaks.iter().map(|(s,x)| (s, x.iter().filter(|y| y.contains(pos)))) {
            for &inner_peak in inner {
                let key = (s, inner_peak, m.pattern_id);
                match pppp.get_mut(&key) {
                    Some(x) => {
                        for haplotype_id in m.haplotype_ids.iter() {
                            x[haplotype_id.sample_id as usize] = x[haplotype_id.sample_id as usize] + 1;
                        }
                    }
                    None => {
                        let mut x = null_count.clone();
                        for haplotype_id in m.haplotype_ids.iter() {
                            x[haplotype_id.sample_id as usize] = x[haplotype_id.sample_id as usize] + 1;
                        }
                        pppp.insert(key, x);
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
    use haplotype::patch_haplotype;

    fn ref_haplotype() -> Vec<NucleotidePos> {
        vec![
            NucleotidePos { nuc: Nucleotide::A, pos: 0 },
            NucleotidePos { nuc: Nucleotide::C, pos: 1 },
            NucleotidePos { nuc: Nucleotide::G, pos: 2 },
            NucleotidePos { nuc: Nucleotide::T, pos: 3 }
        ]
    }

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn test_matches() {
        let c = Weight::new(0, 1000, 0, 0);
        let g = Weight::new(0, 0, 1000, 0);
        let pwm = PWM {weights: vec![c,g], name: "pwm".to_string(), pattern_id: 5, min_score: 1500, direction: PWMDirection::P};
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

        let diffs2 = vec![Diff { pos: 2, reference: nucs("C"), alternative: nucs("NN") }];
        let patched2 = patch_haplotype(&Range::new(1,2), &diffs2, &ref_haplotype());
        let expected2 = vec![nucp('C',1), nucp('N',2), nucp('N',2)];
        assert_eq!(patched2, expected2);

        let diffs3 = vec![Diff { pos: 3, reference: nucs("C"), alternative: nucs("NN") }];
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

    #[test]
    fn test_match_gataa() {
        let pwm = PWM { weights: vec![Weight::new(0,0,100,0), Weight::new(100,0,0,0), Weight::new(0,0,0,100), Weight::new(100,0,0,0), Weight::new(100,0,0,0), ], name: "Example".to_string(), pattern_id: 123, min_score: 499, direction: PWMDirection::P };
        let haplotype_with_padding = vec![nucp('N',0), nucp('G',1), nucp('A',2), nucp('T',3), nucp('A',4), nucp('A',5), nucp('N',6)];
        let haplotype_without_padding = vec![nucp('G',1), nucp('A',2), nucp('T',3), nucp('A',4), nucp('A',5)];
        let haplotype_ids = Rc::new(vec![HaplotypeId { sample_id: 456, side: HaplotypeSide::Right }]);
        let ms = matches(&pwm, &haplotype_with_padding, haplotype_ids.clone());
        assert_eq!(ms.len(), 1);
        let ms2 = matches(&pwm, &haplotype_without_padding, haplotype_ids.clone());
        assert_eq!(ms2.len(), 1);

        let pwm2 = PWM { weights: vec![Weight::new(0,0,100,0), Weight::new(100,0,0,0), Weight::new(0,0,0,100), Weight::new(100,0,0,0), Weight::new(100,0,0,0), ], name: "Example".to_string(), pattern_id: 123, min_score: 500, direction: PWMDirection::P };
        let ms3 = matches(&pwm2, &haplotype_with_padding, haplotype_ids.clone());
        assert_eq!(ms3.len(), 0);
        let ms4 = matches(&pwm2, &haplotype_without_padding, haplotype_ids.clone());
        assert_eq!(ms4.len(), 0);
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