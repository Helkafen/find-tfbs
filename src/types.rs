use std::rc::Rc;
use std::fmt;
use super::range::*;

#[derive(Eq, PartialEq, Clone, Ord, PartialOrd, Debug, Hash, Copy)]
pub enum Nucleotide {
    A, C, G, T, N
}

impl fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Nucleotide::A => write!(f, "A"),
            Nucleotide::C => write!(f, "C"),
            Nucleotide::G => write!(f, "G"),
            Nucleotide::T => write!(f, "T"),
            Nucleotide::N => write!(f, "N"),
        }
    }
}


#[derive(Eq, PartialEq, Ord, PartialOrd, Debug, Clone, Hash)]
pub struct NucleotidePos {
    pub nuc: Nucleotide,
    pub pos: u64
}

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub enum HaplotypeSide { Left, Right }

#[derive(Eq, PartialEq, Debug, Clone)]
pub struct Match {
    pub range: Range,
    pub pattern_id: u16,
    pub haplotype_ids: Rc<Vec<HaplotypeId>>
}

#[derive(Eq, PartialEq, Clone, Ord, PartialOrd, Hash, Debug)]
pub struct Diff {
    pub pos: u64,
    pub reference: Vec<Nucleotide>,
    pub alternative: Vec<Nucleotide>,
}

impl fmt::Display for Diff {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let r = {
            let mut s = "".to_string();
            for x in &self.reference {
                s = s + &format!("{}", x).to_string();
            }
            s
        };
        let a = {
            let mut s = "".to_string();
            for x in &self.alternative {
                s = s + &format!("{}", x).to_string();
            }
            s
        };
        write!(f, "{} {}->{}", self.pos, r, a)
    }
}

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub struct HaplotypeId {
    pub sample_id: usize,
    pub side: HaplotypeSide
}

#[derive(Eq, Clone, PartialEq, Debug)]
pub enum PWMDirection {
    P, N
}

impl fmt::Display for PWMDirection {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            PWMDirection::P => write!(f, "+"),
            PWMDirection::N => write!(f, "-"),
        }
    }
}

#[derive(Eq, PartialEq, Clone, Debug)]
pub enum Pattern {
     PWM { weights: Vec<Weight>, name: String, pattern_id: u16, min_score: i32, direction: PWMDirection },
     OtherPattern {  name: String, pattern_id: u16 }
}

pub fn pattern_length(p: Pattern) -> u32 {
    match p {
        Pattern::PWM{weights, name: _, pattern_id: _, min_score: _, direction: _} => {
            weights.len() as u32
        }
        Pattern::OtherPattern{name: _, pattern_id: _} => {
            0
        }
    }
}

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Weight {
    pub acgtn : Vec<i32>
}

impl Weight {
    pub fn new(a: i32, c: i32, g: i32, t: i32) -> Weight {
        let mut acgtn = vec![a,c,g,t,0];
        acgtn.shrink_to_fit();
        return Weight { acgtn : acgtn };
    }
}