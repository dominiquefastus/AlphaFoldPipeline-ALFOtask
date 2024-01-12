#!/usr/bin/env python3

"""
Script for plotting the pLDDT scores and B-factors for multiple PDB files. 
This script will plot the pLDDT scores and B-factors for each residue in the PDB file
based on a boxplot comparison. It runs like this:

usage: plddt_bval.py [-h] -p ALPHAFOLD_PDB REFERENCE_PDB [-l LABEL]

Process AlphaFold and reference PDB files.

optional arguments:
  -h, --help            show this help message and exit
  -p ALPHAFOLD_PDB REFERENCE_PDB, --pair ALPHAFOLD_PDB REFERENCE_PDB
                        Specify a pair of PDB files
  -l LABEL, --label LABEL
                        Label for each pair
                        
Author:     D. Fastus
"""

import argparse
import sys
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import PDB

def extract_plddt(file_path):
    """
    Extracts pLDDT scores from an AlphaFold PDB file.
    This function should be adapted based on the specific format of your PDB files.
    """
    
    # Placeholder for pLDDT extraction logic
    plddt_scores = np.random.uniform(50, 100, 100)  # Replace with actual extraction logic
    return plddt_scores

def extract_b_factors(file_path):
    """
    Extracts B-factors from a standard PDB file using Biopython.
    """
    
    parser = PDB.PDBParser()
    structure = parser.get_structure('structure', file_path)

    # Extract B-factors from each atom
    # This can be adapted based on the specific format of your PDB files
    b_factors = [atom.get_bfactor() for model in structure for chain in model for residue in chain for atom in residue]

    return b_factors

def generate_box_plot(data, total_pairs, labels=None):
    """
    Generates a single box plot for all pairs, with specific colors, adjusted spacing, and a legend.
    """
    
    # Colors
    colors = [(167/255, 85/255, 104/255), (71/255, 148/255, 149/255)]

    fig, ax = plt.subplots(figsize=(8, 4), dpi=300)
    
    # Calculate positions for each box plot
    # Increase the gap between the pairs
    positions = []
    gap = 2  # Gap between pairs
    for i in range(total_pairs):
        positions.append(1 + (gap + 1) * i)
        positions.append(2 + (gap + 1) * i)

    # Creating box plots with slimmer lines
    box = ax.boxplot(data, patch_artist=False, positions=positions, widths=0.6, showcaps=True, whiskerprops={'linewidth':1.5}, medianprops={'linewidth':1.5}, boxprops={'linewidth':1.5})

    # Setting colors and labels for the legend
    for i, line in enumerate(box['boxes']):
        line.set_color(colors[i % len(colors)])
    legend_labels = ['pLDDT Scores', 'B-Factors']
    ax.legend(handles=[plt.Line2D([0], [0], color=color, lw=4) for color in colors], labels=legend_labels)

    # Customizing plot appearance
    plt.ylabel('Scores/Values')
    
    # Use custom labels if provided, otherwise generate automatic labels
    if labels:
        pair_labels = labels
    else:
        pair_labels = [f'Pair {i + 1}' for i in range(total_pairs)]
    plt.xticks([(1.5 + (gap + 1) * i) for i in range(total_pairs)], pair_labels)

    plt.grid(True)
    plt.savefig('combined_boxplot.png')
    plt.close()

def process_pairs(pairs, labels):
    """
    Process each pair of files and generate a single combined plot.
    """
    
    all_data = []
    # Extract data from each pair
    # generate box plot for each pair
    for af_file, ref_file in pairs:
        af_data = extract_plddt(af_file)
        ref_data = extract_b_factors(ref_file)
        min_length = min(len(af_data), len(ref_data))
        all_data.append(af_data[:min_length])
        all_data.append(ref_data[:min_length])

    generate_box_plot(all_data, len(pairs), labels)

def main():
    # add simple arguments for parser
    parser = argparse.ArgumentParser(description="Process AlphaFold and reference PDB files.")
    parser.add_argument('-p', '--pair', action='append', nargs=2, metavar=('ALPHAFOLD_PDB', 'REFERENCE_PDB'), help='Specify a pair of PDB files', required=True)
    parser.add_argument('-l', '--label', action='append', help='Label for each pair', default=[])

    args = parser.parse_args()

    # Check if the number of labels matches the number of pairs
    if args.label and len(args.label) != len(args.pair):
        parser.error("Number of labels must match the number of pairs")

    # Process each specified pair
    process_pairs(args.pair, args.label)

if __name__ == "__main__":
    main()