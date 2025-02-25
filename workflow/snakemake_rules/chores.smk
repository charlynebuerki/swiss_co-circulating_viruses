'''
This part of the workflow assembles the dataset and corrects 
any discrepencies between local and all sequences
'''

##############################
# Update strain names
###############################

'''rule update_strain_names:
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
'''



##############################
# Add additional sequences
# if you have sequences that are not on NCBI Virus
###############################


## align local sequences:
rule align_local: 
    message:
        """
        Aligning local sequences to {input.reference} using Nextalign
        """
    input:
        sequences = rules.fetch.output.local_sequences,
        reference = rules.copy_reference.output.destination
    output:
        alignment = "data/{strain}/aligned_local.fasta"

    params:
        nuc_mismatch_all = config['align']['allowed_mismatches'],
        nuc_seed_length = 30
    shell:
        """
        nextclade run \
        {input.sequences}  \
        --input-ref {input.reference}\
        --allowed-mismatches {params.nuc_mismatch_all} \
        --min-length {params.nuc_seed_length} \
        --include-reference false \
        --output-fasta {output.alignment} 
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
            --metadata-id-columns strain \
            --output-metadata {output.metadata}
        """

rule add_local_sequences:
    message: 
        """
            merging all aligned sequences into one file. 
        """
    input:
        local_sequences = rules.align_local.output.alignment,
        all_sequences = "data/{strain}/aligned.fasta"
    output:
        sequences = "data/{strain}/all_aligned_sequences.fasta"
    shell:
        """
        cat {input.local_sequences} {input.all_sequences} > {output.sequences}
        """

rule rename_metadata_columns:
    message:
        """
        renaming metadata columns 
        """
    input: 
        meta = rules.add_local_metadata.output.metadata
    output:
        meta_updated = "data/{strain}/all_metadata_updated.tsv"

    params:
        qc_name = "accession",
        qc_old_name = "genbank_accession"

    shell: 
        """
        awk -F'\t' 'NR==1{{gsub("{params.qc_old_name}", "{params.qc_name}"); print; next}} 1' {input.meta} > {output.meta_updated}


        """
