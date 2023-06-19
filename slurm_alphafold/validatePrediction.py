"""
AlphaFold visualization for preassesment

Author:     D. Fastus
"""

# 3d visualization with plddt score gradient

from Bio.PDB import PDBParser, MMCIFParser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import pickle
import json
import requests
import os
import pyrama

class Mode():
    pass

class PDBvisualizer():

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file

    def visualize(self, protein_name, image_path, json_path):
        parser = PDBParser()

        structure = parser.get_structure(f"{protein_name}", self.pdb_file)

# adapted from 
class predicitonScore():

    def plotPLDDT(self, path):
        with open(path,'r') as rank_file:
            ranking = json.load(rank_file)

        model1 = pickle.load(open("/data/staff/biomax/domfas/pipeline/alf_output/1408159/7QRZ/result_model_1.pkl",'rb'))

        # Set the style
        # sns.set(style="darkgrid")

        # Increase figure size and DPI for better resolution
        plt.figure(figsize=(14, 6), dpi=600)

        # Set a custom color for the line
        line_color = '#1f77b4'

        # Plot the data
        plt.plot(model1["plddt"], color=line_color, linewidth=1.5, label="model1")

        # Set plot title and adjust font size
        plt.title("Predicted LDDT per position", fontsize=16)

        # Set axis labels and adjust font size
        plt.xlabel("Positions", fontsize=12)
        plt.ylabel("Predicted LDDT", fontsize=12)

        # Set axis limits and ticks
        plt.ylim(0, 100)
        plt.xlim(0, len(model1["plddt"]))

        # Customize tick parameters
        plt.tick_params(axis='both', which='both', direction='in', bottom=True, top=True, left=True, right=True)

        # Add a grid to the plot
        plt.grid(linestyle='dashed', alpha=0.5)

        # Move the legend to the upper left corner
        plt.legend(loc="lower left", fontsize=12)

        # Save the plot to a file with higher resolution
        plt.savefig("PLDDT.png", dpi=600)

        # Show the plot
        plt.show()
            

    def plotPAE():
        pass

    def plotDISTOGRAM():
        pass

    def plotRamachandran(self, pdb_file=None, pdb_id=None, save_path=None):
        # plot ramachandran plot
        # fetch structure if not provided

        if pdb_file is None and pdb_id is not None:
            pdb = f"https://files.rcsb.org/download/{pdb_id}.pdb"

            response = requests.get(pdb)
            if response.status_code == 200:
                with open(save_path, 'wb') as file:
                    file.write(response.content)
                print(f"PDB file {pdb_id}.pdb downloaded successfully!")
            else:
                print(f"Failed to download PDB file {pdb_id}.pdb")

        elif pdb_file is not None and pdb_id is None:
            pdb = pdb_file
        else:
            raise ValueError("Either pdb_file or pdb_id must be provided")
            sys.exit(1)

        if save_path is None:
            save_path = os.getcwd()

        # plot ramachandran plot with matplotlib using pdb file
        normals, outliers = pyrama.calc_ramachandran(pdb_file)
        pyrama.plot_ramachandran(normals, outliers)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize and validate predicted protein structure')
    parser.add_argument('-p', '--pdb', type=str, help='Path to PDB file')
    parser.add_argument('-id', '--pdb_id', type=str, help='PDB ID')
    parser.add_argument('-s', '--save_path', type=str, help='Path to save plots')
    parser.add_argument('-j', '--json', type=str, help='Path to json file with plddt scores')
    args = parser.parse_args()

    pdb_file = args.pdb
    pdb_id = args.pdb_id
    save_path = args.save_path
    json_path = args.json

    predicitonScore = predicitonScore()
    ramachandran_plot = predicitonScore.plotRamachandran(pdb_file=pdb_file)

