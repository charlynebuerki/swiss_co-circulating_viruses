"""
This part of the workflow handles fetching sequences from various sources.
Uses `config.sources` to determine which sequences to include in final output.

Currently only fetches sequences from GenBank, but other sources can be
defined in the config. If adding other sources, add a new rule upstream
of rule `fetch_all_sequences` to create the file `data/{source}.ndjson` or the
file must exist as a static file in the repo.

Produces final output as

    sequences_ndjson = "data/{strain}/sequences.ndjson"

"""


rule fetch_ncbi_dataset_package:
    params:
        ncbi_taxon_id=lambda w: config["ncbi_taxon_id"].get(w.strain, None),
    output:
        dataset_package=temp("data/{strain}/ncbi_dataset.zip"),
    retries: 5  # Requires snakemake 7.7.0 or later
    benchmark:
        "benchmarks/{strain}/fetch_ncbi_dataset_package.txt"
    shell:
        """
        datasets download virus genome taxon {params.ncbi_taxon_id:q} \
            --no-progressbar \
            --filename {output.dataset_package}
        """


rule extract_ncbi_dataset_sequences:
    input:
        dataset_package="data/{strain}/ncbi_dataset.zip",
    output:
        ncbi_dataset_sequences=temp("data/{strain}/ncbi_dataset_sequences.fasta"),
    benchmark:
        "benchmarks/{strain}/extract_ncbi_dataset_sequences.txt"
    shell:
        """
        unzip -jp {input.dataset_package} \
            ncbi_dataset/data/genomic.fna > {output.ncbi_dataset_sequences}
        """


def _get_ncbi_dataset_field_mnemonics(wildcards) -> str:
    """
    Return list of NCBI Dataset report field mnemonics for fields that we want
    to parse out of the dataset report. The column names in the output TSV
    are different from the mnemonics.

    See NCBI Dataset docs for full list of available fields and their column
    names in the output:
    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/dataformat/tsv/dataformat_tsv_virus-genome/#fields
    """
    fields = [
        "accession",
        "sourcedb",
        "isolate-lineage",
        "geo-region",
        "geo-location",
        "isolate-collection-date",
        "release-date",
        "update-date",
        "length",
        "host-name",
        "isolate-lineage-source",
        "bioprojects",
        "biosample-acc",
        "sra-accs",
        "submitter-names",
        "submitter-affiliation",
    ]
    return ",".join(fields)


rule format_ncbi_dataset_report:
    # Formats the headers to be the same as before we used NCBI Datasets
    # The only fields we do not have equivalents for are "title" and "publications"
    input:
        dataset_package="data/{strain}/ncbi_dataset.zip",
        ncbi_field_map=config["ncbi_field_map"]
    output:
        ncbi_dataset_tsv=temp("data/{strain}/ncbi_dataset_report.tsv"),
    params:
        fields_to_include= _get_ncbi_dataset_field_mnemonics,
    benchmark: "benchmarks/{strain}/format_ncbi_dataset_report.txt"
    shell:
        """
        dataformat tsv virus-genome \
            --package {input.dataset_package} \
            --fields {params.fields_to_include:q} \
            | csvtk -tl rename2 -F -f '*' -p '(.+)' -r '{{kv}}' -k {input.ncbi_field_map} \
            | csvtk -tl mutate -f genbank_accession_rev -n genbank_accession -p "^(.+?)\." \
            | tsv-select -H -f genbank_accession --rest last \
            > {output.ncbi_dataset_tsv}
        """


rule format_ncbi_datasets_ndjson:
    input:
        ncbi_dataset_sequences="data/{strain}/ncbi_dataset_sequences.fasta",
        ncbi_dataset_tsv="data/{strain}/ncbi_dataset_report.tsv",
    output:
        ndjson="data/{strain}/ncbi.ndjson",
    log: "logs/{strain}/format_ncbi_datasets_ndjson.txt",
    benchmark: "benchmarks/{strain}/format_ncbi_datasets_ndjson.txt"
    shell:
        """
        augur curate passthru \
            --metadata {input.ncbi_dataset_tsv} \
            --fasta {input.ncbi_dataset_sequences} \
            --seq-id-column genbank_accession_rev \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            2> {log} > {output.ndjson}
        """

rule fetch_all_genbank_sequences:
    input:
        #all_taxon_ids=expand(
        #    "data/{strain}/ncbi.ndjson", strain=list(config["ncbi_taxon_id"].keys())
        #),
        all_taxon_ids = "data/{strain}/ncbi.ndjson"
    output:
        sequences_ndjson="data/{strain}/genbank.ndjson",
    shell:
        """
        cat {input.all_taxon_ids} > {output.sequences_ndjson}
        """

def _get_all_sources(wildcards):
    return [f"data/{wildcards.strain}/{source}.ndjson" for source in config["sources"]]


rule fetch_all_sequences:
    input:
        all_sources=_get_all_sources,
    output:
        sequences_ndjson="data/{strain}/sequences.ndjson",
    shell:
        """
        cat {input.all_sources} > {output.sequences_ndjson}
        """
