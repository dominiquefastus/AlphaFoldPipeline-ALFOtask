#!/usr/bin/env python3

"""
Script for plotting the B-value per residue for multiple PDB files with optional alignment.

This script will plot the average B-factor for each residue in the PDB file. And runs like this:
usage: bfactor_plot.py [-h] [-l LABELS [LABELS ...]] [-a] pdb_files [pdb_files ...]

Plot B-factors for multiple PDB files.

positional arguments:
  pdb_files             Paths to PDB files

optional arguments:
  -h, --help            show this help message and exit
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
                        Labels for the PDB files
  -a, --align           Align residues across different PDB files

Author:     D. Fastus
"""

import argparse
from biopandas.pdb import PandasPdb
import matplotlib.pyplot as plt
from matplotlib import style
import pandas as pd


def load_and_prepare_pdb(path, align_residues):
    """
    Loads a PDB file and prepare the B-factor data.
    Returns a DataFrame with 'residue_number' and 'b_factor'.
    """
    try:
        # Load the PDB file and extract the B-factors
        pdb = PandasPdb().read_pdb(path)
        df = pdb.df['ATOM'][['residue_number', 'b_factor']]
        grouped = df.groupby('residue_number').mean().reset_index()  # Average B-factor for each residue

        # if the residue numbers don't align, then align them when using the --align flag
        # this is done by subtracting the minimum residue number from all residue numbers
        if align_residues:
            grouped['residue_number'] -= grouped['residue_number'].min()  # Align residue numbering

        return grouped
    
    except Exception as e:
        print(f"Error loading {path}: {e}")
        return pd.DataFrame()


def plot_b_factors(pdb_files, labels, align_residues):
    """
    Plots the B-factors for each PDB file.
    """
    
    plt.figure(figsize=(8, 4), dpi=300)

    # Define the colors for the plots
    colors = [(167/255, 85/255, 104/255), (71/255, 148/255, 149/255)]

    # Plot the B-factors for each PDB file
    for path, label, color in zip(pdb_files, labels, colors):
        pdb_data = load_and_prepare_pdb(path, align_residues)
        if not pdb_data.empty:
            plt.plot(pdb_data['residue_number'], pdb_data['b_factor'], label=label, color=color)

    plt.xlabel('Aligned Residue Number' if align_residues else 'Residue Number')
    plt.ylabel('Average B-Factor')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.savefig('b_factors.png')


def main():
    # add arguments for parser
    parser = argparse.ArgumentParser(description="Plot B-factors for multiple PDB files.")
    parser.add_argument('pdb_files', nargs='+', help='Paths to PDB files')
    parser.add_argument('-l', '--labels', nargs='+', help='Labels for the PDB files', default=[])
    parser.add_argument('-a', '--align', action='store_true', help='Align residues across different PDB files')

    args = parser.parse_args()

    # Check if number of labels and pdb_files match
    if args.labels and len(args.labels) != len(args.pdb_files):
        parser.error("Number of labels must match the number of PDB files.")

    # Use filenames as labels if no labels are provided
    if not args.labels:
        args.labels = [f"File {i+1}" for i in range(len(args.pdb_files))]

    plot_b_factors(args.pdb_files, args.labels, args.align)


if __name__ == "__main__":
    main()
