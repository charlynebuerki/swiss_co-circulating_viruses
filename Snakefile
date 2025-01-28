###############
# Snakemake execution templates:

# To run a default whole genome run ( <6400bp):
# snakemake whole_genome/auspice/cva16_whole_genome.json --cores 1

###############
configfile: "config/configfile.yaml"

wildcard_constraints:
    build="whole",
    strain= r"coronavirus_229E|coronavirus_HKU1|coronavirus_NL63|coronavirus_OC43"
    #gene="|-ncp|-pp|-D1|-D2|-C|-MP|-FP|-HN|-L"
   
#     #from: https://bitbucket.org/snakemake/snakemake/issues/910/empty-wildcard-assignment-works-only-if


build_dir = 'results'
auspice_dir = 'auspice'

# Expand augur JSON paths
rule all:
    input:
        tree_nextclade=expand("results/{strain}/nextclade_dataset_{build}/tree.json", 
                                strain=config.get("strains", ['coronavirus_229E']),
                                build= config.get("builds_to_run", ['whole'])
                                ),
        pathogen_nextclade=expand("results/{strain}/nextclade_dataset_{build}/pathogen.json",
                                strain=config.get("strains", ['coronavirus_229E']),
                                build= config.get("builds_to_run", ['whole'])
                                ),
        sequences_nextclade =expand("results/{strain}/nextclade_dataset_{build}/annotation.gff3",
                                strain=config.get("strains", ['coronavirus_229E']),
                                build= config.get("builds_to_run", ['whole'])
                                ),
        reference_nextclade = expand("results/{strain}/nextclade_dataset_{build}/reference.fasta", 
                                strain=config.get("strains", ['coronavirus_229E']),
                                build= config.get("builds_to_run", ['whole'])
                                )


# Rule to handle configuration files
rule files:
    input:
        #sequence_length =   "{seg}",
        exclude =   config['exclude'],
        lat_longs = config['files']['latitude_longitude_schemes'],
        colors =    config['files']['color_schemes'],
        reference = lambda wildcards: config['files']['reference_config'][wildcards.strain],
        auspice_config = lambda wildcards: config['files']['auspice_config'][wildcards.strain] ,
        clades =   lambda wildcards:      config['files']['clades_config'][wildcards.strain],
        meta=               "data/{strain}/metadata_updated.tsv",
        last_updated_file = "data/{strain}/date_last_updated.txt", #no need
        local_accn_file =   "data/{strain}/local_accn.txt" #no need

files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

include: "workflow/snakemake_rules/import_from_ingest.smk"
include: "workflow/snakemake_rules/chores.smk"
include: "workflow/snakemake_rules/core.smk"


'''
rule update_sequences:
    input:
        sequences = "data/sequences.fasta",
        metadata=files.meta #,
        #add_metadata = files.extended_metafile
    output:
        sequences = "data/sequences_added.fasta"
    params:
        file_ending = "data/all_consensus_hpiv3/*.fa*",
        temp = "data/temp_sequences_added.fasta",
        date_last_updated = files.last_updated_file,
        local_accn = files.local_accn_file, # --dates {params.date_last_updated} --local_accession {params.local_accn} --meta {input.metadata} --add {input.add_metadata}
    shell:
        """
        touch {params.temp} && rm {params.temp}
        cat {params.file_ending} > {params.temp}
        python scripts/update_sequences.py --in_seq {params.temp} --out_seq {output.sequences} 
        rm {params.temp}
        awk '/^>/{{if (seen[$1]++ == 0) print; next}} !/^>/{{print}}' {output.sequences} > {params.temp} && mv {params.temp} {output.sequences}
        """
'''
    
##############################
# Merge all metadata files (NCBI download and own files) and clean them up
# potentially use augur merge: but not the same output can be achieved with augur
###############################
'''
rule add_metadata:
    message:
        """
        Cleaning data in metadata
        """
    input:
        metadata=files.meta,
        new_data=rules.curate_meta_dates.output.final_metadata,
        renamed_strains="data/updated_strain_names.tsv"
    params:
        strain_id_field="accession",
        last_updated = files.last_updated_file,
        local_accn = files.local_accn_file,
    output:
        metadata="data/final_metadata.tsv"
    shell:
        """
        python scripts/add_metadata.py \
            --input {input.metadata} \
            --add {input.new_data} \
            --rename {input.renamed_strains} \
            --local {params.local_accn} \
            --update {params.last_updated}\
            --id {params.strain_id_field} \
            --output {output.metadata}
        """
'''


##############################
# Rest of the augur pipeline
###############################

        

# ##############################
rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice", 
        "data"
    shell:
        "rm -rfv {params}"
