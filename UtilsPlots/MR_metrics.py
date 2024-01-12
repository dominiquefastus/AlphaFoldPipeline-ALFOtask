#!/usr/bin/env python3

"""
Script for plotting the molecular replacement metrics from dimple.log file.
Generally it will plot the R-factor and R-free values over refinement cycles. As well as the TFZ and LLG values over iterations.

The script runs like this:
usage: MR_metrics.py [-h] logfile

Plot data from dimple.log file.

positional arguments:
  logfile     Path to the dimple.log file

optional arguments:
  -h, --help  show this help message and exit

Author:     D. Fastus
"""

import argparse
import re
from configparser import ConfigParser
import matplotlib.pyplot as plt
import numpy as np

def parse_dimple_log_metrics(log_contents):
    """
    Parses the dimple.log file for the B-factor, solvent percentage, R-factor, and R-free values.
    # Function to setup the parser for the metrics in the dimple.log file.
    """
    # Parse the dimple.log file
    # configure the the parser based on the log file
    config = ConfigParser()
    config.read_string(log_contents)

    # parse the B-factor and solvent percentage from the truncate section
    # sometimes the log file does not contain this section like truncate or rwcontents
    b_factor = None
    solvent_percent = None
    if 'truncate' in config:
        try:
            # get the b-factor from the truncate section
            b_factor = float(config['truncate'].get('B-factor', None))
        except Exception:
            print('Could not parse B-factor from truncate section.')
            pass
    if 'rwcontents' in config:
        try:
            # get the solvent percentage from the rwcontents section    
            solvent_percent = float(config['rwcontents'].get('solvent_percent', None))
        except Exception:
            print('Could not parse solvent_percent from rwcontents section.')
            pass

    # parse the R-factor and R-free values from the refmac5 section
    
    # setup dictionary to store the R-factor and R-free values for each refinement type and cycle
    # additional the rms values are stored
    r_values = {'jelly': {'overall_r': [], 'free_r': [], 'rms_bond': [], 'rms_angle': [], 'rms_chiral': []},
                'restr': {'overall_r': [], 'free_r': [], 'rms_bond': [], 'rms_angle': [], 'rms_chiral': []}}

    # parse the R-factor and R-free values from the refmac5 section
    # loop through the refmac5 section and parse the R-factor and R-free values for each refinement type and cycle
    for key in ['refmac5 jelly', 'refmac5 restr']:
        if key in config:
            section = config[key]
            r_values[key.split()[1]] = {
                'overall_r': [float(x) for x in section.get('iter_overall_r', '').strip('[]').split(',')],
                'free_r': [float(x) for x in section.get('iter_free_r', '').strip('[]').split(',')],
                'rms_bond': [float(x) for x in section.get('rmsBOND', '').strip('[]').split(',')],
                'rms_angle': [float(x) for x in section.get('rmsANGL', '').strip('[]').split(',')],
                'rms_chiral': [float(x) for x in section.get('rmsCHIRAL', '').strip('[]').split(',')]
            }

    return b_factor, solvent_percent, r_values

# Function to parse the dimple.log file for the TFZ and LLG values
def parse_dimple_log_phaser(log_contents):
    """
    Parses the dimple.log file for the TFZ and LLG values.
    """
    
    tfz_pattern = r"TFZ==?(\d+\.?\d*)"
    llg_pattern = r"LLG=(\d+)"

    # use regex to parse the TFZ and LLG values from the phaser section
    tfz_values = [float(value) for value in re.findall(tfz_pattern, log_contents)]
    llg_values = [int(value) for value in re.findall(llg_pattern, log_contents)]

    return tfz_values, llg_values

def create_plots(b_factor, solvent_percent, r_values, tfz_values, llg_values):
    """
    Creates the plots from the parsed data.
    """
    
    # B-factor and Solvent Percentage Plot
    # not used in the current version of the script
    '''
    # B-factor and Solvent Percentage Plot
    if b_factor is not None and solvent_percent is not None:
        plt.figure(figsize=(8, 6))
        plt.bar(['B-factor', 'Solvent %'], [b_factor, solvent_percent], color=['skyblue', 'lightgreen'])
        plt.title('B-factor and Solvent Percentage')
        plt.ylabel('Value')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.savefig('b_factor_and_solvent_percentage.png')
    '''

    # R-factor and R-free Values Over Refinement Cycles
    # get the R-factor and R-free values from the parsed data
    # arange the cycles for the x-axis
    jelly_r = r_values['jelly']
    restr_r = r_values['restr']
    overall_r = jelly_r['overall_r'] + restr_r['overall_r']
    free_r = jelly_r['free_r'] + restr_r['free_r']
    cycles = np.arange(1, len(overall_r) + 1)
    
    # plot the R-factor and R-free values over refinement cycles
    plt.figure(figsize=(6, 4), dpi=300)
    plt.plot(cycles, overall_r, label='Overall R-factor', marker='o', color='lightseagreen')
    plt.plot(cycles, free_r, label='Free R-factor', marker='x', color='salmon')

    # plt.title('R-factor and R-free Values Over Refinement Cycles')
    plt.xlabel('Refinement Cycle')
    plt.ylabel('R-factor Value')
    plt.xticks(cycles)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend()
    plt.savefig('r_factors_over_cycles.png')

    # plot the RMS Bond, Angle, and Chiral Values Over Refinement Cycles
    # not used in the current version of the script
    '''
    # RMS Bond, Angle, and Chiral Values Over Refinement Cycles
    rms_bond = jelly_r['rms_bond'] + restr_r['rms_bond']
    rms_angle = jelly_r['rms_angle'] + restr_r['rms_angle']
    rms_chiral = jelly_r['rms_chiral'] + restr_r['rms_chiral']
    plt.figure(figsize=(8, 5))
    plt.plot(cycles, rms_bond, label='RMS Bond Lengths', marker='o', color='mediumseagreen')
    plt.plot(cycles, rms_angle, label='RMS Angles', marker='x', color='mediumpurple')
    plt.plot(cycles, rms_chiral, label='RMS Chiral Centers', marker='s', color='tomato')
    # plt.title('RMS Values Over Refinement Cycles')
    plt.xlabel('Refinement Cycle')
    plt.ylabel('RMS Value')
    plt.xticks(cycles)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.legend()
    plt.savefig('rms_values_over_cycles.png')
    '''

    # TFZ and LLG Values Over Iterations
    # plot if the TFZ and LLG values are available
    if tfz_values and llg_values:
        fig, ax1 = plt.subplots(figsize=(6, 4), dpi=300)
        color = 'seagreen'
        ax1.set_xlabel('Iteration')
        ax1.set_ylabel('TFZ', color=color)
        ax1.plot(tfz_values, color=color, marker='o', label='TFZ')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(axis='y', linestyle='--', alpha=0.7)

        ax2 = ax1.twinx()
        color = 'indianred'
        ax2.set_ylabel('LLG', color=color)
        ax2.plot(llg_values, color=color, marker='x', label='LLG')
        ax2.tick_params(axis='y', labelcolor=color)

        fig.tight_layout()
        # plt.title('TFZ and LLG Values Over Iterations')
        plt.savefig('tfz_llg_values.png')

def main():
    # setup an easy to use command line interface
    parser = argparse.ArgumentParser(description='Plot data from dimple.log file.')
    parser.add_argument('logfile', type=str, help='Path to the dimple.log file')
    args = parser.parse_args()

    with open(args.logfile, 'r') as file:
        log_contents = file.read()

    # retrieve the metrics data from the dimple.log file
    b_factor, solvent_percent, r_values = parse_dimple_log_metrics(log_contents)
    tfz_values, llg_values = parse_dimple_log_phaser(log_contents)

    # create the plots from the parsed data
    create_plots(b_factor, solvent_percent, r_values, tfz_values, llg_values)

if __name__ == "__main__":
    main()