###############
# Snakemake execution templates:

# To run a default whole genome run ( <6400bp):
# snakemake whole_genome/auspice/cva16_whole_genome.json --cores 1

###############
wildcard_constraints:
    seg="whole",
    #gene="|-5utr|-vp4|-vp2|-vp3|-vp1|-2A|-2B|-2C|-3A|-3B|-3C|-3D|-3utr"
   
#     #from: https://bitbucket.org/snakemake/snakemake/issues/910/empty-wildcard-assignment-works-only-if

# Define segments to analyze
segments = ['whole']
#GENES=["-5utr","-vp4", "-vp2", "-vp3", "-vp1", "-2A", "-2B", "-2C", "-3A", "-3B", "-3C", "-3D","-3utr"]
'''GENES = [
    "nucleocapsid_protein",
    "phosphoprotein",
    "D_protein",
    "C_protein",
    "matrix_protein",
    "fusion_protein",
    "hemagglutinin-neuraminidase",
    "large_protein",
]'''
# Expand augur JSON paths
rule all:
    input:
        #augur_jsons = expand("auspice/hpiv3_{segs}.json", segs=segments)
        #augur_jsons = "results/nextclade_dataset/tree.json", "results/nextclade_dataset/sequences"
        tree_nextcalde=expand("results/nextclade_dataset_{segs}/tree.json", segs=segments),
        pathogen_nextclade=expand("results/nextclade_dataset_{segs}/pathogen.json",segs=segments),
        sequences_nextclade =expand("results/nextclade_dataset_{segs}/annotation.gff3",segs=segments),
        reference_nextclade = expand("results/nextclade_dataset_{segs}/reference.fasta", segs=segments)

#rule all_genes:
#    input:
#        augur_jsons = expand("auspice/cva16_whole_genome{genes}.json", genes=GENES)

# Rule to copy over reference genomes from ingest/data/references
rule copy_reference:
    input:
        source="ingest/data/references/reference.fasta"
    output:
        destination="pathogen/config/reference.fasta"
    shell:
        """
        cp {input.source} {output.destination}
        """


# Rule to handle configuration files
rule files:
    input:
        sequence_length =   "{seg}",
        dropped_strains =   "config/dropped_strains.txt",
        lat_longs =         "config/lat_longs.tsv",
        colors =            "config/colors.tsv",
        reference =         "pathogen/config/reference.gb",
        auspice_config =    "pathogen/config/auspice_config.json",
        clades =            "pathogen/clades_genome.tsv",
        meta=               "data/metadata_updated.tsv",
        #extended_metafile=  "data/assign_publications_corrected.tsv",
        last_updated_file = "data/date_last_updated.txt",
        local_accn_file =   "data/local_accn.txt"

files = rules.files.input

##############################
# Download from NBCI Virus with ingest snakefile
###############################

rule fetch:
    input:
        dir = "ingest"
    output:
        sequences="data/sequences.fasta",
        metadata= "data/metadata.tsv",
        local_sequences = "data/local_sequences.fasta",
        local_metadata = "data/local_metadata.fasta"
    params:
        seq="ingest/results/sequences.fasta",
        meta="ingest/results/metadata.tsv",
        seq_loc = "ingest/results/local_sequences.fasta",
        meta_loc = "ingest/data/pathogen_local/1_HPIV_3_meta.tsv"
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
        file_out = "data/updated_strain_names.tsv"
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
        sequences = "data/all_sequences.fasta"
    shell:
        """
        cat {input.local_sequences} {input.all_sequences} > {output.sequences}
        """



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
# Change the format of the dates in the metadata
# Attention: ```augur curate``` only accepts iso 8 formats; please make sure that you save e.g. Excel files in the correct format
###############################

'''rule curate_meta_dates:
    message:
        """
        Cleaning up metadata with augur curate
        """
    input:
        metadata=files.extended_metafile,  # Path to input metadata file
        genbank_meta="data/metadata/genbank_metadata_additional.tsv"  # Generated with bin/extract_genbank_metadata.py
    params:
        strain_id_field="accession",
        date_column="date",
        format=['%Y', '%m.%Y', '%d.%m.%Y', "%b-%Y", "%d-%b-%Y","%Y-%m-%d"],
        temp_metadata="data/temp_curated.tsv"  # Temporary file
    output:
        metadata="data/assign_publications_curated.tsv",  # Final output file for metadata
        genbank_meta="data/metadata/genbank_metadata_curated.tsv",  # Curated genbank metadata
        final_metadata="data/assign_publications_corr_fetched.tsv"  # Final merged output file
    shell:
        """
        # Normalize strings for metadata
        augur curate normalize-strings --metadata {input.metadata} \
            --id-column {params.strain_id_field} \
            --output-metadata {params.temp_metadata}

        # Format dates for metadata
        augur curate format-dates \
            --metadata {params.temp_metadata} \
            --date-fields {params.date_column} \
            --no-mask-failure \
            --expected-date-formats {params.format} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.metadata}
        
        # Remove temporary file
        rm {params.temp_metadata}

        # Normalize strings for genbank metadata
        augur curate normalize-strings --metadata {input.genbank_meta} \
            --id-column {params.strain_id_field} \
            --output-metadata {params.temp_metadata}

        # Format dates for genbank metadata
        augur curate format-dates \
            --metadata {params.temp_metadata} \
            --date-fields {params.date_column} \
            --no-mask-failure \
            --expected-date-formats {params.format} \
            --id-column {params.strain_id_field} \
            --output-metadata {output.genbank_meta}
        
        # Remove temporary file
        rm {params.temp_metadata}

        # Merge curated metadata
        augur merge --metadata meta={output.metadata} genbank_meta={output.genbank_meta}\
            --metadata-id-columns {params.strain_id_field} \
            --output-metadata {output.final_metadata}
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
rule add_local_metadata:
    message:
        """
            adding in local metadata
        """
    input:
        metadata= rules.fetch.output.metadata,
        local_metadata = rules.fetch.output.local_metadata
    output:
        metadata = "data/all_metadata.tsv"
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
        meta_updated = "data/metadata_updated.tsv"

    params:
        qc_name = "qc_overallStatus",
        qc_old_name = "qc.overallStatus"

    shell: 
        """
        awk -F'\t' 'NR==1{{gsub("{params.qc_old_name}", "{params.qc_name}"); print; next}} 1' {input.meta} > {output.meta_updated}

        head -n 1 data/metadata_updated.tsv


        """


##############################
# Rest of the augur pipeline
###############################

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = rules.add_local_sequences.output.sequences
    output:
        sequence_index = "results/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.add_local_sequences.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = files.meta,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        group_by = "country",
        sequences_per_group = 100, # set lower if you want to have a max sequences per group
        strain_id_field= "accession",
        min_date = 1800, 
        min_length = 13000 #length of 15462
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --sequence-index {input.sequence_index} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --exclude {input.exclude} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length} \
            --query 'qc_overallStatus== "good" | qc_overallStatus == "mediocre" | database== "ReVSeq" ' \
            --output {output.sequences}
        """

#useless step I think if ran in ingest initially
rule reference_gb_to_fasta:
    message:
        """
        Converting reference sequence from genbank to fasta format
        """
    input:
        reference = files.reference

    output:
        reference = "results/reference_sequence.fasta"
    run:
        from Bio import SeqIO 
        SeqIO.convert(input.reference, "genbank", output.reference, "fasta")

rule align: 
    message:
        """
        Aligning sequences to {input.reference} using Nextalign
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.reference_gb_to_fasta.output.reference
    output:
        alignment = "results/aligned.fasta"

    params:
        nuc_mismatch_all = 10,
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


rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        # alignment = rules.fix_align_codon.output.alignment,
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    threads: 9
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --nthreads {threads} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree by rerooting and resolving polytomies
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        # alignment = rules.fix_align_codon.output.alignment,
        alignment = rules.align.output.alignment,
        metadata =  files.meta,
    output:
        # tree = "{seg}/results/tree.nwk",
        # node_data = "{seg}/results/branch_lengths.json"
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 3, # was 3
        strain_id_field ="accession",
        # clock_rate = 0.004, # remove for estimation
        # clock_std_dev = 0.0015
        # clock_rate_string = lambda wildcards: f"--clock-rate 0.004 --clock-std-dev 0.0015" if wildcards.gene else ""
        clock_rate_string = "--clock-rate 0.004 --clock-std-dev 0.0015"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            {params.clock_rate_string} \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """
        

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        # alignment = rules.fix_align_codon.output.alignment,
        alignment = rules.align.output.alignment
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --keep-ambiguous\
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts.json"
        # node_data = "{seg}/results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """

rule clades: 
    message: "Assigning clades according to nucleotide mutations"
    input:
        tree=rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = files.clades
    output:
        # clade_data = "{seg}/results/clades.json"
        clade_data = "results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """
rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = files.meta
    output:
        # node_data = "{seg}/results/traits.json"
        node_data = "results/traits.json",
    params:
        traits = "country",
        strain_id_field= "accession"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-node-data {output.node_data} \
            --columns {params.traits} \
            --confidence
        """

rule export:
    message: "Creating auspice JSONs"
    input:
        tree = rules.refine.output.tree,
        metadata = files.meta,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        clades = rules.clades.output.clade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    params:
        strain_id_field= "accession"
    output:
        auspice_json = "auspice/hpiv3_{segs}.json"
        # auspice_json="auspice/cva16_{seg}-accession.json"
        #input {input.clades}
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

rule assemble_dataset:
    input:
        tree= rules.export.output.auspice_json,
        sequences= rules.add_local_sequences.output.sequences
    params:
        directory = "results/nextclade_dataset_{segs}"
    output:
        tree="results/nextclade_dataset_{segs}/tree.json",
        sequences="results/nextclade_dataset_{segs}/sequences.fasta"
    shell:
        """
        mkdir -p {params.directory}
        sed -E 's/\\"node_attrs\\": \\{{/\\"node_attrs\\": \\{{\\"clade_membership\\": \\{{\\"value\\": \\"\\"\\}}, /g' {input.tree} > {output.tree}
        cp {input.sequences} {output.sequences}
        """



rule assemble_nextclade_dataset:
    input:
        pathogen = "ingest/data/references/pathogen.json",
        sequence_annotation = "ingest/data/references/annotation.gff3",
        reference = "ingest/data/references/reference.fasta"
    params:
        directory = rules.assemble_dataset.params.directory
    output:
        pathogen="results/nextclade_dataset_{segs}/pathogen.json",
        sequence_annotation="results/nextclade_dataset_{segs}/annotation.gff3",
        reference="results/nextclade_dataset_{segs}/reference.fasta"
    shell:
        """
        mkdir -p {params.directory}
        cp -u {input.pathogen} {output.pathogen}
        cp -u {input.sequence_annotation} {output.sequence_annotation}
        cp -u {input.reference} {output.reference}
        """

# ##############################
rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice", 
        "data"
    shell:
        "rm -rfv {params}"
