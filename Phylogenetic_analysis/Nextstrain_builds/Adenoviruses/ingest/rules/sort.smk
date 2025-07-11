"""
This part of the workflow handles sorting downloaded sequences and metadata
by aligning them to reference sequences.

It produces output files as

    metadata = "data/{strain}/metadata.tsv"
    sequences = "data/{strain}/sequences.fasta"

"""

rule sort:
    input:
        sequences = rules.curate.output.sequences
    output:
        sequences = "results/{strain}/sequences.fasta"
    shell:
        '''
        seqkit rmdup {input.sequences} > {output}
        '''


rule metadata:
    input:
        metadata = rules.subset_metadata.output.subset_metadata,
        sequences = "results/{strain}/sequences.fasta"
    output:
        metadata = "data/{strain}/metadata_raw.tsv"
    run:
        import pandas as pd
        from Bio import SeqIO

        strains = [s.id for s in SeqIO.parse(input.sequences, 'fasta')]
        d = pd.read_csv(input.metadata, sep='\t', index_col='accession').loc[strains].drop_duplicates()
        d.to_csv(output.metadata, sep='\t')

rule nextclade:
    input:
        sequences = "results/{strain}/sequences.fasta",
        ref = "data/{strain}/reference.fasta"
    output:
        nextclade = "results/{strain}/nextclade.tsv"
    params:
        dataset = "data/{strain}/",
        output_columns = "seqName clade qc.overallScore qc.overallStatus alignmentScore  alignmentStart  alignmentEnd  coverage dynamic"
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
        nextclade = "results/{strain}/nextclade.tsv",
        metadata = "data/{strain}/metadata_raw.tsv"
    output:
        metadata = "results/{strain}/metadata.tsv"
    shell:
        """
        python3 bin/extend-metadata.py --metadata {input.metadata} \
                                       --id-field accession \
                                       --virus-type {wildcards.strain} \
                                       --nextclade {input.nextclade} \
                                       --output {output.metadata}
        """

