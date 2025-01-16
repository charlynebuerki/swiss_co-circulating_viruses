# Ingest pipeline (adapted) for HPIV-3 sequences

This is the adapted ingest pipeline for mpox virus sequences to HPIV-3 sequences.

Of note, one needs the reference genome in data/references/. You also need the annotated gff3 file of that genome. To create it, you can run this script [generate_from_genbank.py](bin/generate_from_genbank.py) manually. Once you have acquired these files, you may run the ingest pipeline. 

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.

## Usage

> NOTE: All command examples assume you are within the `ingest` directory.
> If running commands from the outer `hpiv-3` directory, please replace the `.` with `ingest`

Fetch sequences with

```sh
nextstrain build . data/sequences.ndjson
```

Run the complete ingest pipeline with

```sh
nextstrain build .
```

This will produce three files (within the `ingest` directory):

- `results/metadata.tsv`
- `results/sequences.fasta`
- `results/local_sequences.fasta`

Run the complete ingest pipeline and upload results to AWS S3 with

```sh
nextstrain build . --configfiles build-configs/nextstrain-automation/config.yaml
```

### Including local sequences not from GenBank

By default, this workflow includes local sequences that are not on GenBank. If you wish to remove this part, you can comment out the Snakemake target in the _get_all_targets rule (line number 27) and the 
rule in the include statement (line 61). 

If you wish to keep these, you must include your fasta sequences in `data/pathogen_local`. The pipeline expects the following:
- a metadata file saved as a `.tsv` file. This file must at least have a column called accession which identifies the accession or reference of the sequence.
- all the different fasta files for your sequences. By default, this pipeline reads all the fasta files in `pathogen_local/` folder, will filter the potentially multi-fastsa files on the reference
accession name, and will concatenate all the fastas in one multifasta file in `results/local_sequences.fasta`. 


## Configuration

Configuration takes place in `defaults/config.yaml` by default.
Optional configs for uploading files and Slack notifications are in `build-configs/nextstrain-automation/config.yaml`.

If you wish to run another pathogen, head there to alter the pathogen taxon id number and accession number. 

### Environment Variables

The complete ingest pipeline with AWS S3 uploads and Slack notifications uses the following environment variables:

#### Required

- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`
- `SLACK_TOKEN`
- `SLACK_CHANNELS`

#### Optional

These are optional environment variables used in our automated pipeline for providing detailed Slack notifications.

- `GITHUB_RUN_ID` - provided via [`github.run_id` in a GitHub Action workflow](https://docs.github.com/en/actions/learn-github-actions/contexts#github-context)
- `AWS_BATCH_JOB_ID` - provided via [AWS Batch Job environment variables](https://docs.aws.amazon.com/batch/latest/userguide/job_env_vars.html)

## Input data

### GenBank data

GenBank sequences and metadata are fetched via [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

## `ingest/vendored`

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo) to manage copies of ingest scripts in [ingest/vendored](./vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest).

See [vendored/README.md](vendored/README.md#vendoring) for instructions on how to update
the vendored scripts.
