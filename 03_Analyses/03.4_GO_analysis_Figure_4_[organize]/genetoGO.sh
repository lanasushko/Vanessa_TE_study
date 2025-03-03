#!/bin/bash

# Input and output file paths
input_file=$1
output_file=$2

# Process input file and write to output file
while IFS=$'\t' read -r gene_name go_terms_str; do
    # Split the GO terms by commas
    IFS=',' read -ra go_terms <<< "$go_terms_str"
    
    # Write each gene name and GO term pair to the output file
    for go_term in "${go_terms[@]}"; do
        echo -e "${gene_name}\t${go_term}" >> "$output_file"
    done
done < "$input_file"
