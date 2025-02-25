import sys
import pandas as pd
import requests
import json


if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python script.py <prefix_file> <source_directory_base> <destination_directory>")
        sys.exit(1)

    input_strains = sys.argv[1]
    output_seq = sys.argv[2]
    output_meta = sys.argv[3]

    # Read strains from formatted list
    with open(input_strains, "r") as f:
        strains = json.loads(f.read())  # Convert from JSON list format

    # Fetch sequences
    seq_response = requests.post(
        "https://lapis.cov-spectrum.org/open/v2/sample/alignedNucleotideSequences",
        headers={"accept": "application/json", "Content-Type": "application/json"},
        json={"dataFormat": "FASTA", "strain": strains},
    )

    # Save sequences
    with open(output_seq, "w") as f:
        f.write(seq_response.text)

    # Fetch metadata
    meta_response = requests.post(
        "https://lapis.cov-spectrum.org/open/v2/sample/details",
        headers={"accept": "application/json", "Content-Type": "application/json"},
        json={"strain": strains, "dataFormat": "TSV"},
    )

    # Save metadata with column rename
    with open(output_meta + ".tmp", "w") as f:
        f.write(meta_response.text)

    # Read TSV and rename column
    df = pd.read_csv(output_meta + ".tmp", sep="\t")
    df.rename(columns={"genbankAccession": "genbank_accession"}, inplace=True)
    df.to_csv(output_meta, sep="\t", index=False)

    # Cleanup
    import os
    os.remove(output_meta + ".tmp")
