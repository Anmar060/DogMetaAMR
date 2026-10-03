#!/bin/bash
#$ -l h_rt=48:00:00
#$ -cwd
#$ -q gpu
#$ -l gpu=1
#$ -N dorado_4mC
#$ -m be
#$ -M s1887840@ed.ac.uk

. /etc/profile.d/modules.sh
module load cuda/11.0.2
module load anaconda/2024.02
conda activate samtools_ana202506
#you must have samtools downloaded in an environment


TARGET_DIR=/home/s1887840/scratch/pods_b1/pods_filtered/all
OUTPUT_DIR=/home/s1887840/scratch/fastq
DORADO=/exports/cmvm/eddie/eb/groups/nanopore_users/sw/dorado-1.2.0-linux-x64

time $DORADO/bin/dorado basecaller --kit-name SQK-NBD114-96 -b 0 -v --recursive $DORADO/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_4mC_5mC@v3 $TARGET_DIR > $TARGET_DIR/4mC_5mC_methylation.bam
##time $DORADO/bin/dorado basecaller --kit-name SQK-NBD114-96 -b 0 -v --recursive $DORADO/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0_6mA@v3 $TARGET_DIR > $TARGET_DIR/6mA_methylation.bam

##time samtools fastq -t -c 2 -0 $OUTPUT_DIR/6mA_methylation.fastq.gz $TARGET_DIR/6mA_methylation.bam
time samtools fastq -t -c 2 -0 $OUTPUT_DIR/4mC_5mC_methylation.fastq.gz $TARGET_DIR/4mC_5mC_methylation.bam
