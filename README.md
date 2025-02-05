# Nextstrain repository for Human parainfluenzas (HPIV-1,2,3,4)

This repository provides the pipeline to run the Nextstrain analysis of the family of HPIV viruses on whole-genome. 

It is adapted from [coxsackievirus_a16](https://github.com/hodcroftlab/coxsackievirus_a16) and the official [rsv build](https://github.com/nextstrain/rsv).

# Repository organization 

- [`ingest/`](./ingest) - Download data from GenBank, clean and curate it. It also contains the snakemake file for automatic download and references for alignment 
- snakefile : the computational pipleine managed by Snakemake. 
- [`config/`](./config) - configuration files for the nextstrain run and individual pathogens.
- [`scripts/`](./scripts) - python and add hoc files needed for snakemake processing 
- [`workflow/`](./workflow) - containes the snakemake rules for the workflow


# Running the whole work flow:

First, populate the pathogen folders with sequences from genbank:
You must first acquire the reference genomes in ingest/data/{pathogen}/. You also need the annotated gff3 file of that genome. To create it, you can run this script: [generate_from_genbank.py](bin/generate_from_genbank.py) manually (i.e from a terminal launch python3 ingest/bin/generate_from_genbank.py --reference "NC_003461" --output-dir "ingest/data/HPIV_1" ) with arguments (1) the genbank reference of the sequence and (2) the output directory. The outputs should be stored in ingest/data/{pathogen}. You can find the genbank references in ingest/defaults/config.yaml under accession_reference. 


Run the workflow:
- run command 'nextstrain build .' provided you have installed the [nextstrain tool](https://docs.nextstrain.org/en/latest/install.html) 

To view the build:
- run command 'nextstrain view .'

# Expected outputs:

- an auspice/ folder will be created with invidual builds for each pathogen
- a results/ folder with all gathered files that are created during the workflow, organized by pathogen.

# Prerequisites and repository organization

This workflow will build the number of strains builds defined in ingest/defaults/config.yaml.

This workflow is specially customized to insert your local, unpublished dataset hosted under ingest/data/{pathogen_name}/project_strains/. Specifically, you need:
- a metadata file called '1_metadata.csv' in each {pathogen_name}/project_strains/ folder. The metadata file can contain any data you wish, however the metadata must contain 3 columns : 'accession', 'date', and 'database'. These correspond to the sample name, the date sampled and the project data name (i.e. ReVSeq). Feel free to add more metadata, however we currently support the following visualizations:
    - 'country': sampling country
    - 'region': sampling region (i.e. Europe, Asia, ...)
    - 'location': state or cantonal sub-regions. 
- your consensus fasta files of each sample located in each {pathogen_name}/project_strains/ folder. They should contain the sample name in their file name and end with"consensus.fa" (i.e. Sample_1_consensus.fa)

Additional configurations include:
- changing the name of the virus you are making this build for in the main auspice config files: these are included in each virus subfolder in config/{strain}/auspice_config.json 
- add additional clade notations for each virus in config/{strain}/clades_genome.tsv
 

