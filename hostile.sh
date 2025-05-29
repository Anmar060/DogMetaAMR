#!/bin/bash

#$ -l h_rt=10:00:00
#$ -cwd
#$ -l h_vmem=80G
#$ -N seqkit
#$ -R y
#$ -m be
#$ -M s1887840@ed.ac.uk

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate hostile_012025

TARGET_DIR=/home/s1887840/mazeri_grp/Ana/batch3/demux/ALL/seqkit
OUTPUT_DIR=/home/s1887840/mazeri_grp/Ana/batch3/demux/ALL/clean

for i in barcode04 barcode05 barcode06 barcode07 barcode09 \
barcode10 barcode11 barcode12 barcode13 barcode14 barcode15 \
barcode16 barcode17 barcode18 barcode19 barcode20 barcode21 \
barcode22 barcode23 barcode24 barcode25 barcode26 barcode27 \
barcode28 barcode29 barcode30 barcode31 barcode32 barcode33 \
barcode34 barcode35 barcode36 barcode37 barcode38 barcode39 \
barcode40 barcode41 barcode42 barcode43 barcode44 barcode45 \
barcode46 barcode47 barcode48 barcode50 barcode51 barcode52 ; do

hostile clean --index /home/s1887840/mazeri_grp/Ana/tools/genomes/Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel.fa.gz -o $OUTPUT_DIR --fastq>
done
