#!/bin/bash
#conda activate flye_ana202407

TARGET_DIR=/home/s1887840/Data/dogs/filtered/unduped2
OUTPUT_DIR=/home/s1887840/Data/dogs/assmbl


pwd

for i in barcode49; do

mkdir $OUTPUT_DIR/flye_$i/
print $OUTPUT_DIR/flye_$i/
flye --nano-hq $TARGET_DIR/unduped_$i.fastq.gz -t 26 --out-dir $OUTPUT_DIR/flye_wdog_$i/ --keep-haplotypes --meta

done
