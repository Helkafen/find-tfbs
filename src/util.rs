use super::types::*;
use super::range;

pub fn to_nucleotide(l: u8) -> Nucleotide {
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

pub fn to_nucleotides(letters: &Vec<u8>) -> Vec<Nucleotide> {
    return letters.iter().map(|&l| to_nucleotide(l)).collect();
}

pub fn to_nucleotides_pos(letters: &Vec<u8>, r: &range::Range) -> Vec<NucleotidePos> {
    let nucleotides = to_nucleotides(letters);
    let mut res = Vec::with_capacity(nucleotides.len());
    let mut pos = r.start;
    for n in nucleotides {
        res.push(NucleotidePos {nuc: n, pos: pos });
        pos = pos + 1;
    }
    res
}
