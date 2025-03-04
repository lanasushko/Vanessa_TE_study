# This version contains commandline input option

import argparse

# Hardcoded chromosome file name
CHROMOSOME_FILE = '/tmp/global2/ssushko/Vanessa_proj/2_synteny/alignments/chromosome_correspondance.table'

# Load the chromosome rows into a list of sets
def load_valid_pairs(filename):
    valid_pairs = []
    with open(filename, 'r') as f:
        for line in f:
            # Create a set for the current row of chromosomes
            row = set(line.strip().split('\t'))
            valid_pairs.append(row)
    return valid_pairs

# Check if a pair exists in any row of valid pairs
def is_valid_alignment(chrom1, chrom2, valid_pairs):
    for pair_set in valid_pairs:
        if chrom1 in pair_set and chrom2 in pair_set:
            return True
    return False

# Filter the PAF file
def filter_paf(paf_file, valid_pairs, output_file):
    with open(paf_file, 'r') as paf, open(output_file, 'w') as output:
        for line in paf:
            # Extract the chromosome names from the PAF line
            parts = line.strip().split('\t')
            if len(parts) >= 6:  # Ensure there are enough columns
                chrom1, chrom2 = parts[0], parts[5]  # Query and target names
                # Check if the pair is valid
                if is_valid_alignment(chrom1, chrom2, valid_pairs):
                    output.write(line)

# Main function to execute the filtering
def main():
    parser = argparse.ArgumentParser(description="Filter a PAF file based on chromosome row constraints.")
    parser.add_argument("paf_file", help="Input PAF file to filter")
    parser.add_argument("output_file", help="Output file for filtered alignments")

    args = parser.parse_args()

    # Load chromosome rows and filter the PAF file
    valid_pairs = load_valid_pairs(CHROMOSOME_FILE)
    filter_paf(args.paf_file, valid_pairs, args.output_file)

if __name__ == "__main__":
    main()
