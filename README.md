# DogMetaAMR
Scripts and notes for the dog microbiome metagenomic sequencing and AMR identification project

As of May 2025
1.	Raw read pre-processing: 
  a.	Dorado (0.7.2, 0.8.2, 0.9.0, 0.1.0) for basecalling pod5 reads with super high accuracy basecalling mode (dna_r10.4.1_e8.2_400bps_sup@v5.0.0) 
  b.	Nanofilt v2.8.0 to remove Q < 10 (90% basecalling accuracy) and reads below 800bp
  c.	SeqKit (v 2.8.2) to handle fastq’s and check if there are any duplicated reads 
  d.	NanoPlot for run stats
2.	Data processing & assembly
  a.	Host removal using Hostile (v 1.1.0) and the Roslin reference assembly for labrador retriever dogs (Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel.fa.gz) 
  b.	Read assembly into contigs using Flye (v 2.9.4) 
  c.	Assembly polishing using one round of Medaka (v 2.0.0) and the r1041_e82_400bps_sup_variant_v5.0.0. 
3.	Microbial identification
  a.	Kraken2 (v 2.1.3) at confidence 0.5, using the RefSeqv205 database
  b.	Bracken (v 2.8), 100 kmer length (I need to see if I can build a bigger kmer)
4. AMR detection under construction
  a. AMR detection using AMRplusplus - database to be defined
  b. AMR detection using ResFinder (KMA, NCBI and RefSeq databases)
  c. ABRicate on assembled reads

Under construction
1. Microbial identification using the RefSeqv205 database but on Metaphlan v.4
2. Bryan database (https://doi.org/10.1099/mgen.0.000949) but larger kmer lengths
3. HERRO correction - improves? makes things worse? -loose many reads, it did not make sense for fecal_db but may for RefSeq205
4. Rarefaction curves: ~1/3 of samples are maybe sufficiently sequenced... new types of ARG keep appearing por 1GB clean reads

Tools to look at
1. DAStool (binning correction)
2. R- decontam (remove contamination)
4. R- kraken-biom
5. 

Previous
Based on https://doi.org/10.1099/mgen.0.000802 (Novel canine metagenome Cusco et al 2022)
Used MetaBAT and Kraken2/Bracken to replace DIAMOND and MEGAN-LR, not ideal since data depth was not similar in depth for the rest of the analysis

