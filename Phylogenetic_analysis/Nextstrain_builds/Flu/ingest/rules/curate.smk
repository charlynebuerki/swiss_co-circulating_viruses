"""
This part of the workflow handles curating the data into standardized
formats and expects input file

    sequences_ndjson = "data/sequences.ndjson"

This will produce output files as

    metadata = "data/subset_metadata.tsv"
    sequences = "data/sequences.fasta"

Parameters are expected to be defined in `config.curate`.
"""


rule fetch_general_geolocation_rules:
    output:
        general_geolocation_rules="data/general-geolocation-rules.tsv",
    params:
        geolocation_rules_url=config["curate"]["geolocation_rules_url"],
    shell:
        """
        curl {params.geolocation_rules_url} > {output.general_geolocation_rules}
        """


rule concat_geolocation_rules:
    input:
        general_geolocation_rules="data/general-geolocation-rules.tsv",
        local_geolocation_rules=config["curate"]["local_geolocation_rules"],
    output:
        all_geolocation_rules="data/all-geolocation-rules.tsv",
    shell:
        """
        cat {input.general_geolocation_rules} {input.local_geolocation_rules} >> {output.all_geolocation_rules}
        """


def format_field_map(field_map: dict[str, str]) -> str:
    """
    Format dict to `"key1"="value1" "key2"="value2"...` for use in shell commands.
    """
    return " ".join([f'"{key}"="{value}"' for key, value in field_map.items()])


rule curate:
    input:
        metadata_input ="data/{strain}/metadata_download.tsv",
        all_geolocation_rules="data/all-geolocation-rules.tsv",
        annotations=config["curate"]["annotations"],
    output:
        metadata="data/{strain}/curated_metadata.tsv",
    log:
        "logs/{strain}_curate.txt",
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        strain_regex=config["curate"]["strain_regex"],
        strain_backup_fields=config["curate"]["strain_backup_fields"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],

        annotations_id=config["curate"]["annotations_id"],
        id_field=config["curate"]["id_field"]
    shell:
        """
        augur curate rename \
            --metadata {input.metadata_input} \
            --field-map {params.field_map} \
            --output-metadata {output.metadata}.step1.temp 2>> {log} && \
        
        augur curate apply-geolocation-rules \
            --metadata {output.metadata}.step1.temp \
            --geolocation-rules {input.all_geolocation_rules} \
            --output-metadata {output.metadata}.step2.temp 2>> {log} && \
        
        augur curate apply-record-annotations \
                --metadata {output.metadata}.step2.temp \
                --annotations {input.annotations} \
                --id-field {params.annotations_id} \
                --output-metadata {output.metadata}  2>> {log} && \

            rm {output.metadata}.step1.temp {output.metadata}.step2.temp 

        
        """


rule subset_metadata:
    input:
        metadata="data/{strain}/curated_metadata.tsv",
    output:
        subset_metadata="data/{strain}/metadata_fi.tsv",
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    shell:
        """
        tsv-select -H -f {params.metadata_fields} \
            {input.metadata} > {output.subset_metadata}
        """
