# find-tfbs

A tool to identify rare non-coding variants associated with complex human traits using open chromatin maps and phased whole-genome sequences.

## Inputs

- Phased WGS data in the BCF format
- Genomic coordinates of regions of interest (support for multiple region sets)
- PWM of transcription factors from HOCOMOCO

## VCF output sample

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	INDIVIDUAL1	INDIVIDUAL2	INDIVIDUAL3	INDIVIDUAL4
1	1	regions1.bed,GATA1_HUMAN.H11MO.1.A,100-150	.	.	.	.	COUNTS=2,4;freqs=1/0/3	GT:DS	0|0:0.0	1|1:2.0	1|1:2.0	1|1:2.0
```

A variation in TFBS was found in the region chr1:100-150 (input from regions1.bed). Three individuals have 4 TFBS of the HOCOMOCO pattern GATA1_HUMAN.H11MO.1.A, and one individual only has two TFBS.

The VCF can be used for association testing.

## How to build

For a local computer, the following is sufficient:

```console
cargo build --release
ls target/release/find-tfbs
```

## How to run the tests suite

```console
./run_tests.sh
```

## How to build statically (to make find-tfbs portable to different computers)

Install [rustup](https://www.rust-lang.org/tools/install), then run the following:

```console
./build_static.sh
ls target/x86_64-unknown-linux-musl/release/find-tfbs
```