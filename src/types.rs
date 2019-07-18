use std::rc::Rc;

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
pub enum HaplotypeSide { Left, Right }

#[derive(Eq, PartialEq, Debug)]
pub struct Match {
    pub pos: u64,
    pub pattern_id: u16,
    pub haplotype_ids: Rc<Vec<HaplotypeId>>
}

#[derive(Eq, PartialEq, Clone, Ord, PartialOrd, Hash, Debug)]
pub struct Diff {
    pub pos: u64,
    pub reference: Vec<Nucleotide>,
    pub alternative: Vec<Nucleotide>,
}

#[derive(Eq, PartialEq, Hash, Debug, Clone)]
pub struct HaplotypeId {
    pub sample_id: usize,
    pub side: HaplotypeSide
}

#[derive(Clone)]
pub struct PWM { pub weights: Vec<Weight>, pub name: String, pub pattern_id: u16, pub min_score: i32 }

#[derive(Clone)]
pub struct Weight {
    pub w_a : i32,
    pub w_c : i32,
    pub w_g : i32,
    pub w_t : i32
}

impl Weight {
    pub fn new(a: i32, c: i32, g: i32, t: i32) -> Weight {
        return Weight { w_a : a, w_c : c, w_g : g, w_t : t };
    }
}