#!/bin/bash

# Download the HOCOMOCO database
if [ ! -d "thresholds" ]; then
	wget https://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_thresholds_HUMAN_mono.tar.gz
	tar xzf HOCOMOCOv11_full_thresholds_HUMAN_mono.tar.gz
	mv HOCOMOCOv11_full_thresholds_HUMAN_mono.tar.gz thresholds
fi
if [ ! -f "HOCOMOCOv11_full_pwms_HUMAN_mono.txt" ]; then
	wget https://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_pwms_HUMAN_mono.txt
fi

# Run all find-tfbs tests
cargo test

# To get information about test coverage, install tarpaulin ("cargo install cargo-tarpaulin")  and run "cargo tarpaulin -v"
