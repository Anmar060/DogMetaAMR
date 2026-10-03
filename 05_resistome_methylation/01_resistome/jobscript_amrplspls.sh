#!/bin/bash
module load anaconda/2024.02
conda activate AMR++_env
cd /home/s1887840/mazeri_grp/Ana/tools/AMRplps/AMRplusplus

nextflow run main_AMR++.nf -profile local\
 --pipeline resistome\
 --reads "path/to/your/reads/name*.fastq.gz"\
 --output "path/to/output/folder/" --snp Y

