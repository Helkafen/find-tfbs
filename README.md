# find-tfbs

A tool to identify rare non-coding variants associated with complex human traits using open chromatin maps and phased whole-genome sequences.

## Inputs

- Phased WGS data in the BCF format
- Genomic coordinates of regions of interest (support for multiple region sets)
- PWM of transcription factors from [HOCOMOCO](https://hocomoco11.autosome.ru)

## VCF output sample

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	INDIVIDUAL1	INDIVIDUAL2	INDIVIDUAL3	INDIVIDUAL4
1	1	regions1.bed,GATA1_HUMAN.H11MO.1.A,100-150	.	.	.	.	COUNTS=2,4;freqs=1/0/3	GT:DS	0|0:0.0	1|1:2.0	1|1:2.0	1|1:2.0
```

A variation in TFBS was found in the region chr1:100-150 (input from regions1.bed). Three individuals have 4 TFBS of the HOCOMOCO pattern GATA1_HUMAN.H11MO.1.A, and one individual only has two TFBS.

The VCF can be used for association testing.

## How to build

Install [rustup](https://www.rust-lang.org/tools/install), and the clang compiler (using your package manager), then run the following in find-tfbs's root directory:

```console
cargo build --release
ls target/release/find-tfbs
```

## How to run the tests suite

```console
./run_tests.sh
```

## How to build statically (to make find-tfbs portable to different computers)

Install [rustup](https://www.rust-lang.org/tools/install) and clang, then run the following:

```console
./build_static.sh
ls target/x86_64-unknown-linux-musl/release/find-tfbs
```

## Usage

```console
find-tfbs
    --input chr12.bcf \
    --chromosome 12 \
    --bed Erythro.bed,Nkcell.bed \
    --output output.vcf.gz \
    --pwm_file HOCOMOCOv11_full_pwms_HUMAN_mono.txt \
    --pwm_names GATA1_HUMAN.H11MO.1.A,GATA2_HUMAN.H11MO.1.A \
    --pwm_threshold 0.0001 \
    --pwm_threshold_directory thresholds \
    --reference hg38.fa
```


## Adapting find-tfbs to use alternatives to PWMs

The data type Pattern (in src/types.rs) can be a PWM or any other pattern type. You can extend this data type and write the corresponding matching algorithm in src/pattern.rs. You will also need to load pattern definitions from a file.

find-tfbs can scan different kinds of patterns simultaneously. The only condition is that all patterns need to have a unique pattern\_id.
