"""
This part of the workflow handles fetching sequences from local repository. Given that the files are multi-fasta files, they must be curated and filtered. You can remove this part of the 
workflow by removing the line in the  ingest/snakemake that includes this workflow. If doing so, you must also comment out in the snakemake rule "all" to remove the target of local sequences. 

expected input of this part:
- fasta files for each sample in data/local_pathogen

output:
- concatenated fasta file for all of the local sequences in "results/{strain}/local_sequences.fasta". 

"""
import re
import glob
import os


# Generate samples dynamically for each strain

samples = {
    strain: [
        re.sub(r"_consensus$", "", os.path.splitext(os.path.basename(f))[0])
        for f in glob.glob(f"data/{strain}/project_strains/*.fa")
    ]
    for strain in strains
}

# Precompute the expanded file paths for all strain/sample combinations
def expand_filtered_paths():
    return [
        f"results/{strain}/filtered_fasta/{segment}/{sample}_filtered.fasta"
        for strain in strains
        for segment in segments
        for sample in samples[strain]
    ]

def expand_modified_paths(segment):
    return [
        f"results/{strain}/filtered_fasta/{segment}/{sample}_modified.fasta"
        for strain in strains
        for sample in samples[strain]
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
        sequences = "data/{strain}/project_strains/{sample}_consensus.fa"
    output:
        filtered_sequence = "results/{strain}/filtered_fasta/{segment}/{sample}_filtered.fasta"
    params:
        pathogen_accession = lambda wildcards: config["local"]["accession_reference"][wildcards.strain][wildcards.segment]
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
        sequences = "results/{strain}/filtered_fasta/{segment}/{sample}_filtered.fasta"
    output:
        sequences = "results/{strain}/filtered_fasta/{segment}/{sample}_modified.fasta"
    shell:
        """
        python bin/rename_fasta_headers.py --input {input.sequences} --output {output.sequences} --name {wildcards.sample}
        """

rule merge_all_fastas:
    message:
        """
        Combining all fasta files in one file
        """
    input:
        lambda wildcards: expand_modified_paths(wildcards.segment)
    output:
        sequences = "results/{strain}/local_sequences_{segment}.fasta"
    params:
        strain="{strain}"
    shell:
        """
        cat {input} > {output.sequences}
        """