#!/bin/bash

TARGET_DIR=/home/s1887840/Data/dogs/filtered
DBNAME=/home/s1887840/shared/data/Kraken_databases/kraken_standard_050624
OUTPUT_DIR=/home/s1887840/Data/dogs/ident

for i in barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 \
barcode07 barcode09 barcode10 barcode11 barcode12 \
barcode13 barcode14 barcode15 barcode16 barcode17 barcode18 \
barcode19 barcode20 barcode21 barcode22 barcode23 barcode24 \
barcode25 barcode26 barcode27 barcode28 barcode29 barcode30 \
barcode31 barcode32 barcode33 barcode34 barcode35 barcode36 \
barcode37 barcode38 barcode39 barcode40 barcode41 barcode42 \
barcode43 barcode44 barcode45 barcode46 barcode47 barcode48 \
barcode49 barcode50 barcode51 barcode52 barcode53 ; do

#mkdir $OUTPUT_DIR/$i
#kraken2 --threads 4 --db $DBNAME --gzip-compressed $TARGET_DIR/unduped_$i.fastq.gz --confidence 0.5 --output $OUTPUT_DIR
#kraken2 --threads 16 --db /home/s1887840/shared/data/Kraken_databases/kraken_standard_050624 --gzip-compressed $TARGET_>
#echo $OUTPUT_DIR/$i
#gzip -c /home/s1887840/Data/dogs/filtered/unduped_$i.fastq.gz > /home/s1887840/Data/dogs/zip/unduped_$i.fastq.gz
echo /home/s1887840/Data/dogs/AMR/barcode_$i.csv >> AMRpBarcode.txt
grep "=," /home/s1887840/Data/dogs/AMR/barcode_$i.csv | wc -l >> AMRpBarcode.txt
done
