'''
This part of the workflow expects input files
            expects local sequences (project sequences) in :
            - data/project_strains/

            with the following structure:
            project_strains: 
            - metadata fiel called : '1_metadata.csv' containing at least 3 columns: 'accession', 'date', and 'database'. These correspond to the sample name, the date sampled and the project data name (i.e. ReVSeq).
            - for location and country and region fields to work properly, metadata file should contain columns: "country", "region", "location". 
            - consensus fasta files of each sample with the sample name in their file name and end with"consensus.fa" (i.e. Sample_1_consensus.fa)

'''

import re
import glob
import os


# Generate samples dynamically for each strain

samples = {
        re.sub(r"_consensus$", "", os.path.splitext(os.path.basename(f))[0])
        for f in glob.glob(f"data/project_strains/*.fa")
}

# Precompute the expanded file paths for all strain/sample combinations
def expand_filtered_paths():
    return [
        f"data/filtered_fasta/{sample}_filtered.fasta"
        for sample in samples
    ]

def expand_modified_paths():
    return [
        f"data/filtered_fasta/{sample}_modified.fasta"
        for sample in samples
    ]



rule filter_all_fastas:
    input:
        expand_filtered_paths()


rule filter_on_pathogen:
    message: 
        """
        Fetching local sequences and filtering on pathogen accession 
        """
    input:
        sequences = "data/project_strains/{sample}_consensus.fa"
    output:
        filtered_sequence = "data/filtered_fasta/{sample}_filtered.fasta"
    params:
        pathogen_accession = lambda wildcards: config["local"]["accession_reference"]
    shell:
        """
        seqkit grep -r -n -p "{params.pathogen_accession}" {input.sequences} > {output.filtered_sequence}
        """

rule curate_fasta_header:
    message:
        """
        Replacing accession by sample name.
        """
    input:
        sequences = "data/filtered_fasta/{sample}_filtered.fasta"
    output:
        sequences = "data/filtered_fasta/{sample}_modified.fasta"
    shell:
        """
        python scripts/rename_fasta_headers.py --input {input.sequences} --output {output.sequences} --name {wildcards.sample}
        """

rule merge_all_fastas:
    message:
        """
        Combining all fasta files in one file
        """
    input:
        expand_modified_paths()
    output:
        sequences = "data/local_sequences.fasta"
    shell:
        """
        cat {input} > {output.sequences}
        """


rule add_local_sequences:
    input:
        local_sequences = "data/local_sequences.fasta",
        all_sequences = "data/sequences.fasta"
    output:
        sequences = "data/all_sequences.fasta"
    shell:
        """
        cat {input.local_sequences} {input.all_sequences} > {output.sequences}
        """

rule add_local_metadata:
    message:
        """
            adding in local metadata
        """
    input:
        metadata= "data/metadata.tsv",
        local_metadata = "data/project_strains/1_metadata.csv"
    output:
        metadata = "data/all_metadata.tsv"
    shell: 
        """
        augur merge \
            --metadata metadata={input.metadata} \
            --metadata local_metadata={input.local_metadata} \
            --metadata-id-columns accession \
            --output-metadata {output.metadata}
        """
