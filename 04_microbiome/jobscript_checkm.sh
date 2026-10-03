#!/bin/bash

TARGET_DIR=/home/s1887840/Data/dogs/assmbl/medaka
COV=/home/s1887840/Data/dogs/assmbl/medaka
OUTPUT_DIR=/home/s1887840/Data/dogs/metabat2

for i in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 \
barcode07 barcode08 barcode09 barcode10 barcode11 barcode12 \
barcode13 barcode14 barcode15 barcode16 barcode17 barcode18 \
barcode19 barcode20 barcode21 barcode22 barcode23 barcode24 \
barcode25 barcode26 barcode27 barcode28 barcode29 barcode30 \
barcode31 barcode32 barcode33 barcode34 barcode35 barcode36 \
barcode37 barcode38 barcode39 barcode40 barcode41 barcode42 \
barcode43 barcode44 barcode45 barcode46 barcode47 barcode48 \
barcode49 barcode50 barcode51 barcode52 barcode53 ; do

#metabat2 -i $TARGET_DIR/medaka_$i/consensus.fasta -o $OUTPUT_DIR/metabat2_$i/ -v
#checkm lineage_wf -f /home/s1887840/Data/dogs/metabat2/CheckM_barcode14.txt -t 8 -x fa /home/s1887840/Data/dogs/metabat>
checkm lineage_wf -f $OUTPUT_DIR/CheckM_$i.txt -t 8 -x fa $OUTPUT_DIR/metabat2_$i $OUTPUT_DIR/metabat2_$i/SCG_$i
done
