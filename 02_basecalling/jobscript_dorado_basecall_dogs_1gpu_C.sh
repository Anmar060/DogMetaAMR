#!/bin/bash
#$ -l h_rt=48:00:00
#$ -cwd
#$ -q gpu
#$ -pe gpu-a100 1
#$ -l h_vmem=80G
#$ -N dorado_1gpu_C
#$ -R y
#$ -m be
#$ -M s1887840@ed.ac.uk
 
. /etc/profile.d/modules.sh
module load cuda/11.0.2
 
time /home/s1887840/mazeri_grp/Ana/tools/dorado-0.7.1-linux-x64/bin/dorado basecaller --kit-name SQK-NBD114-96 -b 0 -v --recursive /home/s1887840/mazeri_grp/Ana/tools/dorado-0.7.1-linux-x64/bin/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 /home/s1887840/mazeri_grp/Ana/DogMetaAMR_1C > DogMetaAMR_1C_1cpu.bam



