'''
core nextclade / nextstrain workflow
'''

rule index_sequences:
    message:
        """
        Creating an index of sequence composition for filtering
        """
    input:
        sequences = rules.add_local_sequences.output.sequences
    output:
        sequence_index = build_dir + "/{strain}/{build}/sequence_index.tsv"
    shell:
        """
        augur index \
            --sequences {input.sequences} \
            --output {output.sequence_index}
        """

rule filter_recent:
    message:
        """
        Filtering recent sequences to: 
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.add_local_sequences.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = files.meta,
        exclude = files.exclude
    output:
        sequences = build_dir +"/{strain}/{build}/filtered_recent.fasta",
        log = build_dir + "/{strain}/{build}/filtered_recent.log"
    params:
        group_by = config['filter']['group_by']['recent'],
        sequences_per_group = config['filter']['subsample_max_sequences']['recent'], 
        strain_id_field= "accession",
        min_date = lambda wildcards: config['filter']['resolution']['min_date'][wildcards.strain], 
        max_date = lambda wildcards: config['filter']['resolution']['max_date'],
        min_length = lambda wildcards: config['filter']['min_length'][wildcards.strain][wildcards.build],
        min_coverage = f"genome_coverage>{config['filter']['min_coverage']}",
        include = "config/{strain}/strains_to_include.txt"
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
            --max-date {params.max_date} \
            --query '{params.min_coverage} & bioproject_accession != "PRJEB83635" & host =="Homo sapiens" | database== "ReVSeq" ' \
            --include-where 'database="ReVSeq"' \
            --include {params.include} \
            --output {output.sequences} \
            --output-log {output.log}
        """

rule filter_background:
    message:
        """
        Filtering background sequences to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.add_local_sequences.output.sequences,
        sequence_index = rules.index_sequences.output.sequence_index,
        metadata = files.meta,
        exclude = files.exclude
    output:
        sequences = build_dir +"/{strain}/{build}/filtered_background.fasta",
        log = build_dir + "/{strain}/{build}/filtered_background.log"
    params:
        group_by = config['filter']['group_by']['background'],
        sequences_per_group = config['filter']['subsample_max_sequences']['background'], 
        strain_id_field= "accession",
        min_date = lambda wildcards: config['filter']['resolution']['background_min_date'][wildcards.strain], 
        max_date = lambda wildcards: config['filter']['resolution']['min_date'][wildcards.strain],
        min_length = lambda wildcards: config['filter']['min_length'][wildcards.strain][wildcards.build],
        min_coverage = f"genome_coverage>{config['filter']['min_coverage']}"
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
            --max-date {params.max_date} \
            --min-length {params.min_length} \
            --query '{params.min_coverage} & bioproject_accession != "PRJEB83635" & host =="Homo sapiens" ' \
            --output {output.sequences} \
            --output-log {output.log}
        """

rule combine_samples:
    input:
        subsamples = [rules.filter_recent.output.sequences, rules.filter_background.output.sequences] 
    output:
        sequences = build_dir + "/{strain}/{build}/filtered.fasta"
    shell:
        """
        cat {input.subsamples} | seqkit rmdup > {output}
        """

rule align: 
    message:
        """
        Aligning sequences to {input.reference} using Nextalign
        """
    input:
        sequences = rules.combine_samples.output.sequences,
        reference = rules.copy_reference.output.destination
    output:
        alignment = "results/{strain}/{build}/aligned.fasta"

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

rule tree:
    message:
        """
        Creating a maximum likelihood tree
        """
    input:
        # alignment = rules.fix_align_codon.output.alignment,
        alignment = rules.align.output.alignment
    output:
        tree = build_dir+"/{strain}/{build}/tree_raw.nwk"
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
        alignment = rules.align.output.alignment,
        metadata =  files.meta,
    output:
        tree = build_dir+"/{strain}/{build}/tree.nwk",
        node_data = build_dir+"/{strain}/{build}/branch_lengths.json"
    params:
        root = lambda wildcards: config['refine']['root'][wildcards.strain][wildcards.build],
        coalescent = config['refine']['coalescent'],
        date_inference = config['refine']['date_inference'],
        clock_filter_iqd = config['refine']['clock_filter_iqd'],
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
            --root {params.root} \
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
        alignment = rules.align.output.alignment
    output:
        node_data = build_dir+"/{strain}/{build}/nt_muts.json"
    params:
        inference = config['ancestral']['inference']
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
        node_data = build_dir+"/{strain}/{build}/aa_muts.json"
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
        clade_data = build_dir+"/{strain}/{build}/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """

rule subclades: 
    message: "Assigning subclades according to nucleotide mutations"
    input:
        tree=rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        subclades = files.subclades
    params:
        membership_name = "subclade_membership"
    output:
        subclade_data = build_dir+"/{strain}/{build}/subclades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.subclades} \
            --membership-name {params.membership_name} \
            --output-node-data {output.subclade_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.traits!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = files.meta
    output:
        # node_data = "{seg}/results/traits.json"
        node_data = build_dir+"/{strain}/{build}/traits.json",
    params:
        traits = config['traits']['columns'],
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
        subclades = rules.subclades.output.subclade_data,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    params:
        strain_id_field= "accession"
    output:
        auspice_json = auspice_dir+"/{strain}_{build}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} {input.subclades} \
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
        directory = build_dir+"/{strain}/nextclade_dataset_{build}"
    output:
        tree= build_dir+"/{strain}/nextclade_dataset_{build}/tree.json",
        sequences= build_dir+"/{strain}/nextclade_dataset_{build}/sequences.fasta"
    shell:
        """
        mkdir -p {params.directory}
        sed -E 's/\\"node_attrs\\": \\{{/\\"node_attrs\\": \\{{\\"clade_membership\\": \\{{\\"value\\": \\"\\"\\}}, /g' {input.tree} > {output.tree}
        cp {input.sequences} {output.sequences}
        """

rule assemble_nextclade_dataset:
    input:
        pathogen = lambda wildcards: config['files']['pathogen_config'][wildcards.strain][wildcards.build],
        sequence_annotation = lambda wildcards: config['files']['sequence_annotation'][wildcards.strain][wildcards.build],
        reference =  rules.copy_reference.output.destination
    params:
        directory = rules.assemble_dataset.params.directory
    output:
        pathogen=build_dir+"/{strain}/nextclade_dataset_{build}/pathogen.json",
        sequence_annotation=build_dir+"/{strain}/nextclade_dataset_{build}/annotation.gff3",
        reference=build_dir+"/{strain}/nextclade_dataset_{build}/reference.fasta"
    shell:
        """
        mkdir -p {params.directory}
        cp -u {input.pathogen} {output.pathogen}
        cp -u {input.sequence_annotation} {output.sequence_annotation}
        cp -u {input.reference} {output.reference}
        """