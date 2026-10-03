#!/bin/bash
#$ -l h_rt=48:00:00
#$ -cwd
#$ -q gpu
#$ -pe gpu-a100 1
#$ -l h_vmem=80G
#$ -N dorado_demux_D
#$ -R y
#$ -m be
#$ -M s1887840@ed.ac.uk
 
. /etc/profile.d/modules.sh
module load cuda/11.0.2
 
time /home/s1887840/mazeri_grp/Ana/tools/dorado-0.7.1-linux-x64/bin/dorado demux --no-classify -v --output-dir /home/s1887840/mazeri_grp/Ana/demux /home/s1887840/mazeri_grp/Ana/DogMetaAMR_1D_1cpu.bam



