"""
This part of the workflow handles fetching sequences from local repository. Given that the files are multi-fasta files, they must be curated and filtered. You can remove this part of the 
workflow by removing the line in snakemake that includes this workflow. If doing so, you must also comment out in the snakemake rule "all" to remove the target of local sequences. 

expected input of this part:
- fasta files for each sample in data/local_pathogen

output:
- concatenated fasta file for all of the local sequences. 

"""
import glob
import os


#fetch accession reference from the config file 
accession_reference = config["local"]["HPIV_3"]["accession_reference"]

# Dynamically fetch all file names from the specified directory
input_dir = "data/pathogen_local"
samples = [re.sub(r"_consensus$", "", os.path.splitext(os.path.basename(f))[0]) for f in glob.glob(f"{input_dir}/*.fa")]

rule filter_all_fastas:
    input: 
        expand("results/filtered_fasta/{sample}_filtered.fasta", sample=samples)


rule filter_on_pathogen:
    message: 
        """
        fetching local sequences and filtering on pathogen accession 
        """
    input:
        sequences = "data/pathogen_local/{sample}_consensus.fa"
    output:
        filtered_sequence = "results/filtered_fasta/{sample}_filtered.fasta"

    params:
        pathogen_accession = accession_reference
    shell:
        """
        seqkit grep -r -n -p "{params.pathogen_accession}" {input.sequences} > {output.filtered_sequence}

        """

rule curate_fasta_header:
    message:
        """
        replacing accession by sample name.
        """
    input:
        sequences = rules.filter_on_pathogen.output
    output:
        sequences = "results/filtered_fasta/{sample}_modified.fasta"
    shell:
        """
        python bin/rename_fasta_headers.py --input {input.sequences} --output {output.sequences} --name {wildcards.sample}
        """


rule merge_all_fastas:
    message:
        """
        combining all fasta files in one file
        """
    input: 
        expand("results/filtered_fasta/{sample}_modified.fasta", sample=samples)
    output:
        sequences = "results/local_sequences.fasta"
    shell:
        """
        cat {input} > {output.sequences}
        """