#!/usr/bin/env python3

"""
Script for plotting the molecular replacement metrics overview from overall_metrics_results.csv file.
It will plot the TFZ against LLG values for each protein. As well as the clashscore, Ramachandran outliers, R-free, and rotamer outliers.
The gradient plot is similar to the pdb validation plot.

The script runs like this:
usage: plddt_bval.py [-h] -p ALPHAFOLD_PDB REFERENCE_PDB [-l LABEL]

Process AlphaFold and reference PDB files.

optional arguments:
  -h, --help            show this help message and exit
  -p ALPHAFOLD_PDB REFERENCE_PDB, --pair ALPHAFOLD_PDB REFERENCE_PDB
                        Specify a pair of PDB files
  -l LABEL, --label LABEL
                        Label for each pair
(edna2) (base) [domfas@fe1 UtilsPlots]$ python overall_MR.py -h
usage: overall_MR.py [-h] csv_file

Generate protein metrics plots from CSV data.

positional arguments:
  csv_file    Path to the CSV file containing protein metrics.

optional arguments:
  -h, --help  show this help message and exit
                        
Author:     D. Fastus
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from adjustText import adjust_text
import argparse

import math

def read_data(file_path):
    """
    Reads the CSV file and returns the data.
    """
    
    return pd.read_csv(file_path)

def prepare_data(data):
    """
    Transposes and cleans the data for analysis.
    """
    
    # Transpose the data and drop empty columns
    data_transposed = data.transpose()
    data_transposed.columns = data_transposed.iloc[0]
    data_transposed = data_transposed.drop(data_transposed.index[0])
    return data_transposed.dropna(axis=1, how='all')


# Function to create the scatter plot with adjusted label sizes
def create_scatter_plot(data):
    """
    Create and save the scatter plot with larger data point labels.
    """
    
    # Create the scatter plot
    scatter_data = data[['Type', 'LLG', 'TFZ']]
    scatter_data = scatter_data.dropna()
    
    # Convert the LLG and TFZ values to numeric
    scatter_data['LLG'] = pd.to_numeric(scatter_data['LLG'], errors='coerce')
    scatter_data['TFZ'] = pd.to_numeric(scatter_data['TFZ'], errors='coerce')

    # Create the scatter plot
    plt.figure(figsize=(11, 8), dpi=300)
    ax = plt.subplot(111)
    types = scatter_data['Type'].unique()
    
    # loop through the types and plot the data
    # loop through the data and add the protein name as a label
    for t in types:
        subset = scatter_data[scatter_data['Type'] == t]
        plt.scatter(subset['TFZ'], subset['LLG'], s=20, label=t, alpha=0.6)
        for i in subset.index:
            plt.text(subset['TFZ'][i], subset['LLG'][i], i, fontsize=12, ha='right', va='bottom')  

    # setup line to visualize the TFZ and LLG thresholds to solve molecular replacement
    plt.axhline(y=40, color='r', linestyle='-', label='LLG 40 (minimum for correct)')
    plt.axhline(y=60, color='orange', linestyle='-', label='LLG 60 (difficult problems)')
    plt.axhline(y=120, color='green', linestyle='-', label='LLG 120 (ideal minimum)')
    
    plt.axvline(x=5, color='r', linestyle='-', label='TFZ 5 (minimum for solved)')
    plt.axvline(x=6, color='orange', linestyle='-', label='TFZ 6 (possibly solved)')
    plt.axvline(x=8, color='green', linestyle='-', label='TFZ 8 (definitely solved)')
    
    plt.xlabel('TFZ Score', fontsize=14)
    plt.ylabel('LLG Score', fontsize=14)
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fancybox=True, shadow=True, ncol=3,
              fontsize=12)

    plt.grid(True)
    plt.savefig('tfz_vs_llg.png')


def create_gradient_plot(data):
    """
    Create and save the gradient plot.
    """
    
    def draw_custom_gradient_bar_with_corrected_scales(ax, metric_values, y, cmap, label, worst_val, best_val):
        """
        Draws a custom gradient bar with corrected scales.
        """
        
        # create a costum gradient bar with corrected scales for the given metric values
        gradient = np.linspace(0, 1, 256)
        gradient = np.vstack((gradient, gradient))
        ax.imshow(gradient, aspect='auto', cmap=cmap, extent=(0.05, 0.95, y - 0.1, y + 0.1))

        # sort the metric values and define the label positions
        sorted_values = sorted(metric_values.items(), key=lambda x: x[1])
        label_positions = np.linspace(y + 0.28, y + 0.82, len(metric_values))

        # loop through the sorted values and plot the data
        # plots the label and indicator for each protein on the metric gradient bar
        for (protein, value), label_y in zip(sorted_values, label_positions):
            value = 0.95 * value
            ax.plot(value, y, 'v', color='black', markersize=4)
            ax.text(value, label_y - 0.08, f'{protein}', va='center', ha='center', fontsize=14, color='black',
                    bbox=dict(edgecolor='none', fc=(1.,1.,1.), pad=-0.08))
            ax.plot([value, value], [y, label_y - 0.15], color='black', linestyle='-', linewidth=0.6)

        # add descriptions to the gradient bar
        ax.text(0.06, y - 0.2, f'Worst ({worst_val})', va='center', ha='center', fontsize=14, color='black')
        ax.text(0.94, y - 0.2, f'Best ({best_val})', va='center', ha='center', fontsize=14, color='black')
        ax.text(0.5, y + 0.7, label, va='center', ha='center', fontsize=12, color='black', fontweight='bold')

    # set two colors for cusstom gradient bar
    custom_cmap_new = LinearSegmentedColormap.from_list("custom_cmap_new", ["tomato", "mediumseagreen"])

    # create a dictionary with the metrics and their values
    # can be adapted to include more metrics
    metrics = {
        'Clashscore': data['Clashscore'].dropna().astype(float),
        'Ramachandran outliers': data['Ramachandran outliers'].dropna().str.rstrip('%').astype(float),
        'R-free': data['R-free'].dropna().astype(float),
        'Rotamer outliers': data['Rotamer outliers'].dropna().str.rstrip('%').astype(float)
    }

    # normalize the metrics to a scale from 0 to 1 for the gradient plot 
    # the gradients bar have the same length, but the scales are different
    normalized_metrics = {}
    normalized_metrics['R-free'] = ((metrics['R-free'] - 0.6) / (0 - 0.6)).to_dict()
    normalized_metrics['Clashscore'] = ((metrics['Clashscore'] - 20) / (0 - 20)).to_dict()
    normalized_metrics['Ramachandran outliers'] = ((metrics['Ramachandran outliers'] - 3) / (0 - 3)).to_dict()
    normalized_metrics['Rotamer outliers'] = ((metrics['Rotamer outliers'] - 5) / (0 - 5)).to_dict()

    # plot each gradient bar for each metric and viosualize the protein names on the scale
    fig, ax = plt.subplots(figsize=(14, 11), dpi=500)
    draw_custom_gradient_bar_with_corrected_scales(ax, normalized_metrics['R-free'], 4.4, custom_cmap_new, 'R-free', "0.6", "0")
    draw_custom_gradient_bar_with_corrected_scales(ax, normalized_metrics['Clashscore'], 0.8
                                                   , custom_cmap_new, 'Clashscore', "20", "0")
    draw_custom_gradient_bar_with_corrected_scales(ax, normalized_metrics['Ramachandran outliers'], 3.2, custom_cmap_new, 'Ramachandran Outliers', "3%", "0%")
    draw_custom_gradient_bar_with_corrected_scales(ax, normalized_metrics['Rotamer outliers'], 2, custom_cmap_new, 'Rotamer Outliers', "5%", "0%")

    # ax.set_title('Protein Quality Metrics with Custom Gradient Scale')
    ax.set_xticks([0.05, 0.95])
    ax.set_xticklabels(['Worst', 'Best'], fontdict={'fontsize': 14})
    ax.set_xlim(0, 1)
    ax.set_ylim(0.3, 5.4)
    ax.set_yticks([])
    ax.set_xlabel('Metric Scale', fontdict={'fontsize': 14})
    ax.set_ylabel('Metrics', fontdict={'fontsize': 14})

    # add the frame and ticks from the plot
    for spine in ax.spines.values():
        spine.set_visible(True)
        
    
    ax.tick_params(axis=u'both', which=u'both', length=0)
    plt.savefig('custom_gradient_plot.png')


def main():
    # setup an easy to use command line interface
    parser = argparse.ArgumentParser(description="Generate protein metrics plots from CSV data.")
    parser.add_argument("csv_file", help="Path to the CSV file containing protein metrics.")
    args = parser.parse_args()

    data = read_data(args.csv_file)
    prepared_data = prepare_data(data)

    create_scatter_plot(prepared_data)
    create_gradient_plot(prepared_data)
    print("Plots have been generated and saved.")

if __name__ == "__main__":
    main()
