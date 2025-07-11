'''
This part of the workflow assembles the dataset and corrects 
any discrepencies between local and all sequences
'''

##############################
# Update strain names
###############################

rule update_strain_names:
    message:
        """
        Updating strain information in metadata.
        """
    input:
        file_in =  files.meta
    output:
        file_out = "data/{strain}/updated_strain_names.tsv"
    shell:
        """
        time bash scripts/update_strain.sh {input.file_in} {output.file_out}
        """

##############################
# Add additional sequences
# if you have sequences that are not on NCBI Virus
###############################
rule add_local_sequences:
    input:
        local_sequences = rules.fetch.output.local_sequences,
        all_sequences = rules.fetch.output.sequences
    output:
        sequences = "data/{strain}/all_sequences.fasta"
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
        metadata= rules.fetch.output.metadata,
        local_metadata = rules.fetch.output.local_metadata
    output:
        metadata = "data/{strain}/all_metadata.tsv"
    shell: 
        """
        augur merge \
            --metadata metadata={input.metadata} \
            --metadata local_metadata={input.local_metadata} \
            --metadata-id-columns accession \
            --output-metadata {output.metadata}
        """


rule rename_metadata_columns:
    message:
        """
        renaming metadata columns 
        """
    input: 
        meta = rules.add_local_metadata.output.metadata
    output:
        meta_updated = "data/{strain}/metadata_updated.tsv"

    params:
        qc_name = "qc_overallStatus",
        qc_old_name = "qc.overallStatus"

    shell: 
        """
        awk -F'\t' 'NR==1{{gsub("{params.qc_old_name}", "{params.qc_name}"); print; next}} 1' {input.meta} > {output.meta_updated}


        """
