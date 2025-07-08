"""
This part of the workflow handles sorting downloaded sequences and metadata
by aligning them to reference sequences.

It produces output files as

    metadata = "data/{strain}/metadata.tsv"
    sequences = "data/{strain}/sequences.fasta"

"""



#might want to pre-sort the metadata such that only samples with segments of enough length are retained
#would also want to keep only samples with high enough QC from nextclade -> save as columns for each segment in metadata


rule metadata:
    input:
        metadata = rules.subset_metadata.output.subset_metadata,
        sequences = "data/{strain}/sequences_{segment}.fasta"
    output:
        metadata = "data/{strain}/metadata_raw_{segment}.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, 'fasta')]
        d = pd.read_csv(input.metadata, sep='\t', index_col='accession_2').loc[strains].drop_duplicates()
        d.to_csv(output.metadata, sep='\t')

rule nextclade:
    input:
        sequences = "data/{strain}/sequences_{segment}.fasta",
        ref = "data/{strain}/{segment}/reference.fasta" #need to get & change for every build
    output:
        nextclade = "results/{strain}/nextclade_{segment}.tsv"
    params:
        dataset = "data/{strain}/{segment}/",
        output_columns = "seqName clade qc.overallScore qc.overallStatus alignmentScore  alignmentStart  alignmentEnd  coverage dynamic" #need to change such that dynamically it gets a new segment name
    threads: 8
    shell:
        """
        nextclade3 run -D {params.dataset}  -j {threads} \
                          --output-columns-selection {params.output_columns} \
                          --output-tsv {output.nextclade} \
                          {input.sequences}
        """


rule extend_metadata: 
    input:
        nextclade = "results/{strain}/nextclade_{segment}.tsv", 
        metadata = "data/{strain}/metadata_raw_{segment}.tsv"
    output:
        metadata = "results/{strain}/metadata_{segment}.tsv"
    shell:
        """
        python3 bin/extend-metadata.py --metadata {input.metadata} \
                                       --id-field accession_2 \
                                       --virus-type {wildcards.strain} \
                                       --nextclade {input.nextclade} \
                                       --output {output.metadata}
        """

