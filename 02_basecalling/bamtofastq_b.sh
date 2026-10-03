#!/bin/bash

TARGET_DIR=/home/s1887840/mazeri_grp/Ana/demux/1C/demux
OUTPUT_DIR=/home/s1887840/mazeri_grp/Ana/demux/1C/fastq
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

samtools fastq -t --threads 8 -c 2 -0 $OUTPUT_DIR/$i.fastq.gz $TARGET_DIR/SQK-NBD114-96_$i.bam

done
