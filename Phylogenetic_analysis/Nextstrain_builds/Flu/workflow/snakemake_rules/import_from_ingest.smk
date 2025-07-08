'''
This part of the workflow imports the data fetched (local and from NCBI) from ingest & 
ensures correct format. 
it expects input files:
- sequences stored in ingest/results/<pathogen_name>/
- reference.gbk file stored in config/<pathogen_name>/<segment>/
- local sequences stored in ingest/results/<pathogen_name>/
- associated metadata stored in ingest/results/<pathogen_name>/
- local metadata associated with project strains in ingest/<pathogen_name>/project_strains/1_<pathogen_name>_metadata.csv

outputs:
- sequences in data/<pathogen_name>/
- associated metadata in data/<pathogen_name/ (local and non-local metadata)
- local sequences in data/<pathogen_name>/
'''

rule copy_reference:
    input:
        source="ingest/data/{strain}/{build}/reference.fasta"
    output:
        destination="config/{strain}/{build}/reference.fa"
    shell:
        """
        cp {input.source} {output.destination}
        """

rule copy_gff_reference:
    input:
        source="ingest/data/{strain}/{build}/annotation.gff3"
    output:
        destination="config/{strain}/{build}/annotation.gff3"
    shell:
        """
        cp {input.source} {output.destination}
        """


rule fetch:
    input:
        dir = "ingest"
    output:
        sequences="data/{strain}/sequences_{build}.fasta",
        metadata= "data/{strain}/metadata_{build}.tsv",
        local_sequences = "data/{strain}/local_sequences_{build}.fasta"
    params:
        seq="ingest/data/{strain}/sequences_{build}.fasta",
        meta="ingest/results/{strain}/metadata_{build}.tsv",
        seq_loc = "ingest/results/{strain}/local_sequences_{build}.fasta"
    shell:
        """
        cd {input.dir} 
        snakemake --cores 9 all
        cd ../
        cp -u {params.seq} {output.sequences}
        cp -u {params.meta} {output.metadata}
        cp -u {params.seq_loc} {output.local_sequences}
        """