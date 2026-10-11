#!/bin/sh
#$ -l h_rt=48:00:00
#$ -cwd
#$ -pe sharedmem 16
#$ -l h_vmem=16G
#$ -R y
#$ -m be
#$ -M s1887840@ed.ac.uk

. /etc/profile.d/modules.sh
module load anaconda/2024.02
conda activate kraken2_213

INPUT_DIR=/home/s1887840/scratch/merge
OUTPUT_DIR=/home/s1887840/scratch/merge/nanoplots
DBNAME=/home/s1887840/scratch/database/references/kraken_db

for i in 610961 613700 610761 612458 610221 615944 602262 \
609099 618137 615418 618280 610016 616593 612781 \
606949 605133 612994 611670 610585 611399 612225 \
604686 617829 602222 612896 616418 616075 621936 \
602260 617612 604687 617525 608316 615366 604035 \
616954 609569 614260 618278 611235 607416 605007 \
616229 616594 610526 617743 604076 620640 617224 \
604688 606755 605005 603261 608315 615943 610615 \
611871 616290 613581 611579 614991 611234 613128 \
616291 622320 609385 610449 605132 602271 \
605007 616073 614994 610879 607344; do

for l in  01 02 03; do

mkdir $OUTPUT_DIR/"$i"_"$l"/
time kraken2 \
            --threads 16 \
            --confidence 0.5 \
            --report-minimizer-data \
            --db $DBNAME \
            --report $OUTPUT_DIR/"$i"_"$l"/k2_noplant_"$i"_"l".c0.5.report \
            --output $OUTPUT_DIR/"$i"_"$l"/k2_noplant_"$i"_"l".c0.5.kraken \
            --memory-mapping \
            --gzip-compressed \
            $INPUT_DIR/"$i"_"l"*.fastq.gz

time kraken2 \
      --threads 16 \
            --confidence 0.15 \
            --report-minimizer-data \
            --db $DBNAME \
            --report $OUTPUT_DIR/"$i"_"$l"/k2_noplant_"$i"_"l".c0.15.report \
            --output $OUTPUT_DIR/"$i"_"$l"/k2_noplant_"$i"_"l".c0.15.kraken \
            --memory-mapping \
            --gzip-compressed \
            $INPUT_DIR/"$i"_"l"*.fastq.gz

time bracken  -d $DBNAME -i $OUTPUT_DIR/"$i"_"$l"/k2_noplant_"$i"_"$l".c0.5.report -o $OUTPUT_DIR/"$i"_"$l"/brkn_noplant_c05_"$i"_"$l".txt -r 150
time bracken  -d $DBNAME -i $OUTPUT_DIR/"$i"_"$l"/k2_noplant_"$i"_"$l".c0.5.report -o $OUTPUT_DIR/"$i"_"$l"/brkn_noplant_c05_"$i"_"$l".txt -r 15000
time bracken  -d $DBNAME -i $OUTPUT_DIR/"$i"_"$l"/k2_noplant_"$i"_"$l".c0.15.report -o $OUTPUT_DIR/"$i"_"$l"/brkn_noplant_c15_"$i"_"$l".txt -r 150
time bracken  -d $DBNAME -i $OUTPUT_DIR/"$i"_"$l"/k2_noplant_"$i"_"$l".c0.15.report -o $OUTPUT_DIR/"$i"_"$l"/brkn_noplant_c15_"$i"_"$l".txt -r 15000
