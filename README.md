# DogMetaAMR
Research code and computational analysis scripts accompanying my PhD thesis on long-read metagenomic approaches for microbiome, antimicrobial resistance and pathogen detection.

## About this repository

This repository contains the scripts, analysis code and computational notes used across the thesis. *It is intended as a record of the analyses performed and the parameters used, rather than as a standalone production pipeline.* The tools, pipelines nor software belongs to me: the intention of this repository is to make the computational methodology underlying the thesis inspectable and reproducible rather than having large sections of code on my thesis. 

The code was developed iteratively over the course of the project as datasets, research questions and sequencing/software versions changed. Therefore, this repository includes both final analysis scripts and 'intermediary' scripts documenting methodological decisions.

Where possible, final thesis analyses have been consolidated into reproducible scripts with defined inputs and outputs. Some HPC scripts retain University of Edinburgh (Eddie)-specific paths, environments or scheduler settings and would require adaptation before use on another system.



##Repository status: thesis analysis archive, finalised 2026.

###Thesis computational analyses

1. Long-read preprocessing and quality control

Oxford Nanopore sequencing data were processed using combinations of:

Dorado for basecalling and modified-base calling

NanoFilt v2.8.0 for removal of reads with Q < 10 and reads <800 bp

SeqKit v2.8.2 for FASTQ manipulation and duplicate-read checks

NanoPlot / NanoPack for sequencing QC and run statistics

Hostile v1.1.0 for host-read removal using the Canis lupus familiaris ROS_Cfam_1.0 reference genome

Dorado versions changed during the project as ONT software was updated. Final analyses used the basecalling model appropriate to the R10.4.1/E8.2 chemistry and SUP basecalling.


2. Canine faecal metagenome assembly and genome recovery

Long-read metagenomic assemblies were generated using:

Flye v2.9.4 in metagenomic mode (MetaFlye)

one round of Medaka v2.0.0 polishing using r1041_e82_400bps_sup_variant_v5.0.0

MetaBAT2 for metagenomic binning

CheckM for assessment of MAG completeness and contamination

Assembly methods were initially informed by the canine long-read metagenome study of Cuscó et al. (2022):
https://doi.org/10.1099/mgen.0.000802


3. Microbiome taxonomic profiling

Taxonomic classification of cleaned metagenomic reads used:

Kraken2 v2.1.3

confidence threshold 0.5

RefSeq-based Kraken database

Bracken v2.8 for abundance re-estimation

These outputs formed the basis of the longitudinal canine gut microbiome analyses in R.

The final R analysis includes:

taxonomic abundance summaries

alpha diversity

beta diversity / ordination

longitudinal models

diet, antibiotic and gastrointestinal-health associations

differential-abundance analyses

Experimental database and classifier comparisons carried out during method development are retained separately from the final analysis where relevant.


4. Antimicrobial resistance / resistome analysis

The longitudinal resistome analysis used AMR++ / MEGARes, including SNP-confirmed resistance-marker output.

Downstream Chapter 3 analyses include:

sequencing-depth normalisation

ARG abundance and richness

resistance category and subcategory summaries

prevalence and core/accessory resistome structure

longitudinal ARG persistence and turnover

antibiotic, diet and gastrointestinal-health associations

microbiome–resistome relationships

targeted investigation of selected resistance determinants

The final statistical analyses are consolidated into a single Chapter 3 R script.

Assembly-level ARG screening of selected samples was additionally performed using:

ABRicate v1.0.1 (initially for method development, later for corroboration across tools)

MEGARes database

This was used to connect ARG detections to assembled contigs and genome bins rather than as a replacement for the read-level AMR++ analysis.


5. Native-DNA methylation and ARG contextualisation

A subset of four contrasting faecal metagenomes was used for an exploratory analysis of whether native Oxford Nanopore methylation could provide additional ARG or host context.

This analysis included:

Dorado modified-base calling for 6mA and 4mC/5mC

methylation-aware mapping with minimap2

position-level methylation summaries using Modkit v0.6.4

genome-associated motif discovery using NanoMotif v1.1.2

comparison of NanoMotif results with MetaBAT2 bins and CheckM quality

assembly-level ARG contextualisation with ABRicate

Two NanoMotif parameter regimes were retained to assess robustness:

a FLUPOULmeta-inspired parameter set

a more conservative parameter set

The analysis was designed as a feasibility study rather than as a benchmark or optimised methylation pipeline.

A second exploratory read-level analysis, adapted from the CUPID concept, was performed on selected samples. This combined:

NanoMotif motif discovery

Epimetheus read-level methylation profiling

AMRFinderPlus v4.2.7

Kraken2 v2.1.3

MOB-suite v3.1.9

methylation-profile similarity between plasmid-like ARG reads and candidate chromosome reads

These analyses are retained in the repository because they document both the successful and resolution-limited outcomes of methylation-assisted ARG contextualisation.

CUPID, from Ürel (2026) https://doi.org/10.64898/2026.02.16.706128
FLUPOULmeta from Bloemen (2025) https://doi.org/10.3389/fmicb.2025.1614301


6. Statistical analysis

Statistical analyses and visualisation were performed primarily in R.

Final chapter analyses have been consolidated where possible so that each chapter can be reproduced from a single R session rather than from multiple exploratory scripts. This consolidation was carried out using GPT 4.1 on ELM: this was carried out since the smaller scripts are full of personal notes, including dog IDs and personal communications, and branching experiments not represented in the final analysis.

The R code includes combinations of:

linear and generalised linear mixed models

compositional microbiome analyses

PERMANOVA

longitudinal within-animal comparisons

diversity analyses

ARG prevalence and abundance summaries

publication-resolution figures and thesis tables

###Previous
Based on https://doi.org/10.1099/mgen.0.000802 (Novel canine metagenome Cusco et al 2022)
1.) Used MetaBAT and Kraken2/Bracken to replace DIAMOND and MEGAN-LR, not ideal since data depth was not similar in depth for the rest of the analysis.

2.) alternative taxonomic databases

3.) MetaPhlAn comparisons

4.) HERRO read correction

5.) alternative bin-refinement approaches

6.) sequencing-depth / rarefaction investigations

These are retained only where they help document analytical decisions (Appendix) and should not be interpreted as components of the final workflow.

