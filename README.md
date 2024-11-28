# Nextstrain repository for Human parainfluenza 3 virus

This repository provides the pipeline to run the Nextstrain analysis of HPIV-3 on whole-genome. 

It is adapted from coxsackievirus_a16

# Repository organization 

- [`ingest/`](./ingest) - Download data from GenBank, clean and curate it. It also contains the snakemake file for automatic download and references for alignment 
- snakefile : the computational pipleine managed by Snakemake. 
- [`config/`](./config) - configuration files for the nextstrain run general
- [`pathogen/`](./pathogen) -- sequences and configurations files for the whole genome run
- [`scripts/`](./scripts) - python and add hoc files needed for snakemake processing 


make sure that the .gb file is inside of config/pathogen




