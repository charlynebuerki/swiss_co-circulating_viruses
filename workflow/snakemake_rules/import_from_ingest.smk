'''
This part of the workflow imports the data fetched (local and from NCBI) from ingest & 
ensures correct format. 
it expects input files:
- sequences stored in ingest/results/<pathogen_name>/
- reference.gb file stored in ingest/data/<pathogen_name>
- local sequences stored in ingest/results/<pathogen_name>/
- associated metadata stored in ingest/results/<pathogen_name>/
- local metadata associated with project strains in ingest/<pathogen_name>/project_strains/1_<pathogen_name>_metadata.csv

outputs:
- sequences in data/<pathogen_name>/
- associated metadata in data/<pathogen_name/ (local and non-local metadata)
- local sequences in data/<pathogen_name>/
'''


# Rule to copy over reference genomes from ingest/data/references
rule copy_reference:
    input:
        source="ingest/data/{strain}/reference.fasta"
    output:
        destination="config/{strain}/reference.fa"
    shell:
        """
        cp {input.source} {output.destination}
        """

rule copy_genbank_reference:
    input:
        source="ingest/data/{strain}/reference.gbk"
    output:
        destination="config/{strain}/reference.gb"
    shell:
        """
        cp {input.source} {output.destination}
        """

'''rule reference_gb_to_fasta:
    message:
        """
        Converting reference sequence from genbank to fasta format
        """
    input:
        reference = rules.copy_reference.output.destination
    
    params:
        reference_file=lambda wildcards: config['files']['reference_config'][wildcards.strain].split("/")[-1]

    output:
        reference = build_dir+"/{strain}/{params.reference_file}"
    run:
        from Bio import SeqIO 
        SeqIO.convert(input.reference, "genbank", output.reference, "fasta")
'''

rule fetch:
    input:
        dir = "ingest"
    output:
        sequences = "data/{strain}/aligned.fasta",
        metadata= "data/{strain}/metadata.tsv",
        local_sequences = "data/{strain}/local_sequences.fasta",
        local_metadata = "data/{strain}/local_metadata.tsv"
    params:
        seq = "ingest/results/{strain}/aligned.fasta",
        meta="ingest/results/{strain}/metadata.tsv",
        seq_loc = "ingest/results/{strain}/local_sequences.fasta",
        meta_loc = "ingest/data/{strain}/project_strains/1_metadata.tsv"
    shell:
        """
        cd {input.dir} 
        snakemake --cores 9 all
        cd ../
        cp -u {params.seq} {output.sequences}
        cp -u {params.meta} {output.metadata}
        cp -u {params.seq_loc} {output.local_sequences}
        cp -u {params.meta_loc} {output.local_metadata}

        """