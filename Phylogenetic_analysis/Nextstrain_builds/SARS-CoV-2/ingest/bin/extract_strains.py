import sys
import pandas as pd


if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python script.py <prefix_file> <source_directory_base> <destination_directory>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_strains = sys.argv[2]
    output_strains_formatted = sys.argv[3]


    # Read the first column from the TSV file
    df = pd.read_csv(input_file, sep='\t', usecols=[0])

    # Extract strain names
    strains = df.iloc[:, 0].tolist()

    # Save raw strain list
    with open(output_strains, "w") as f:
        f.write("\n".join(strains) + "\n")

    # Format as a JSON-like list
    formatted_strains = "[" + ",\n".join(f'"{strain}"' for strain in strains[1:]) + "]"

    # Save formatted strain list
    with open(output_strains_formatted, "w") as f:
        f.write(formatted_strains)
