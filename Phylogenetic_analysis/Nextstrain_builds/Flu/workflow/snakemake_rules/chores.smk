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

#some housekeeping on the metadata
def format_field_map(field_map: dict[str, str]) -> str:
    """
    Format dict to `"key1"="value1" "key2"="value2"...` for use in shell commands.
    """
    return " ".join([f'"{key}"="{value}"' for key, value in field_map.items()])

rule merge_local_metadata:
    message:
        """
        renaming metadata columns & merging with local metadata
        """
    input: 
        metadata = rules.fetch.output.metadata,
        local_metadata = files.local_meta
    output:
        metadata = "data/{strain}/all_metadata_{build}.tsv"
    log:
        "logs/{strain}_{build}_rename_meta.txt",
    params:
        field_map = format_field_map(config["curate"]["field_map"]),
        field_map_2 = format_field_map(config["curate"]["field_map_2"]),
    shell: 
        """
        augur curate rename \
            --metadata {input.metadata} \
            --field-map {params.field_map} \
            --output-metadata {output.metadata}.step_1.temp 2>> {log} && \
            
        augur curate rename \
            --metadata {output.metadata}.step_1.temp \
            --field-map {params.field_map_2} \
            --output-metadata {output.metadata}.temp 2>> {log} && \

        augur merge \
            --metadata metadata={output.metadata}.temp \
            --metadata local_metadata={input.local_metadata} \
            --metadata-id-columns accession \
            --output-metadata {output.metadata} 2>> {log} && \
        
        rm {output.metadata}.step_1.temp {output.metadata}.temp

        """

rule add_local_sequences:
    input:
        local_sequences = rules.fetch.output.local_sequences,
        all_sequences = rules.fetch.output.sequences
    output:
        sequences = "data/{strain}/all_sequences_{build}.fasta"
    shell:
        """
        cat {input.local_sequences} {input.all_sequences} > {output.sequences}
        """



