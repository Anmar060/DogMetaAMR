#!/bin/bash

TARGET_DIR=/home/s1887840/mazeri_grp/Ana/demux/ALL/filtered
OUTPUT_DIR=/home/s1887840/mazeri_grp/Ana/demux/ALL/dupes

#module load anaconda/2024.02
#conda activate seqkit_ana202407

pwd

#rm /home/s1887840/mazeri_grp/Ana/demux/ALL/dupes/duped_*.fastq.gz

for i in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 \
barcode07 barcode08 barcode09 barcode10 barcode11 barcode12 \
barcode13 barcode14 barcode15 barcode16 barcode17 barcode18 \
barcode19 barcode20 barcode21 barcode22 barcode23 barcode24 \
barcode25 barcode26 barcode27 barcode28 barcode29 barcode30 \
barcode31 barcode32 barcode33 barcode34 barcode35 barcode36 \
barcode37 barcode38 barcode39 barcode40 barcode41 barcode42 \
barcode43 barcode44 barcode45 barcode46 barcode47 barcode48 \
barcode49 barcode50 barcode51 barcode52 barcode53; do


#flye --nano-raw $TARGET_DIR/filtered_$i.fastq.gz --out-dir $OUTPUT_DIR/$i_flye/ --keep-haplotypes --meta --plasmids
seqkit rmdup -s < $TARGET_DIR/filtered_$i.fastq.gz > $OUTPUT_DIR/duped_$i.fastq.gz
done
