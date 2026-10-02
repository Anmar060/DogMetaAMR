#!/bin/bash

TARGET_DIR_AMBL=/home/s1887840/Data/dogs/assmbl
TARGET_DIR_FIL=/home/s1887840/Data/dogs/filtered/unduped2
OUTPUT_DIR=/home/s1887840/Data/dogs/assmbl/medaka

pwd

for i in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 \
barcode07 barcode08 barcode09 barcode10 barcode11 barcode12 \
barcode13 barcode14 barcode15 barcode16 barcode17 barcode18 \
barcode19 barcode20 barcode21 barcode22 barcode23 barcode24 \
barcode25 barcode26 barcode27 barcode28 barcode29 barcode30 \
barcode31 barcode32 barcode33 barcode34 barcode35 barcode36 \
barcode37 barcode38 barcode39 barcode40 barcode41 barcode42 \
barcode43 barcode44 barcode45 barcode46 barcode47 barcode48 \
barcode49 barcode50 barcode51 barcode52 barcode53; do

mkdir $OUTPUT_DIR/medaka_$i/
print $OUTPUT_DIR
#flye --nano-hq $TARGET_DIR/unduped_$i.fastq.gz -t 20 --out-dir $OUTPUT_DIR/flye_$i/ --keep-haplotypes --meta
medaka_consensus -t 20 -i $TARGET_DIR_FIL/unduped_$i.fastq.gz -d $TARGET_DIR_AMBL/flye_$i/assembly.fasta -o $OUTPUT_DIR/>
done 