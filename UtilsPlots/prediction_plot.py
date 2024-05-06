""" 
This is the script to produce AlphaFold prediction quality metrics plots.

The general code for the plots was adapted and modified from:
https://github.com/jasperzuallaert/VIBFold/blob/main/visualize_alphafold_results.py

The adapted code is also featured in:
https://elearning.vib.be/courses/alphafold/ 
https://blog.biostrand.be/explained-how-to-plot-the-prediction-quality-metrics-with-alphafold2

"""

import glob
import math
import os
import numpy as np
from matplotlib import pyplot as plt
import argparse
import pickle

def get_pae_plddt(model_names, is_multimer):
    out = {}
    for i,name in enumerate(model_names):
        d = pickle.load(open(name,'rb'))
        basename = os.path.basename(name)
        basename = basename[basename.index('model'):]
        if is_multimer:
            out[f'{basename}'] = {'plddt': d['plddt'], 'pae':d['predicted_aligned_error']}
        else:
            out[f'{basename}'] = {'plddt': d['plddt']}
    return out

def generate_output_images(feature_dict, out_dir, name, pae_plddt_per_model, is_multimer):
    msa = feature_dict['msa']
    seqid = (np.array(msa[0] == msa).mean(-1))
    seqid_sort = seqid.argsort()
    non_gaps = (msa != 21).astype(float)
    non_gaps[non_gaps == 0] = np.nan
    final = non_gaps[seqid_sort] * seqid[seqid_sort, None]

    ##################################################################
    plt.figure(figsize=(9, 6), dpi=300)
    ##################################################################
    plt.title(f"Sequence coverage ({name})")
    plt.imshow(final,
               interpolation='nearest', aspect='auto',
               cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
    plt.plot((msa != 21).sum(0), color='black')
    plt.xlim(-0.5, msa.shape[1] - 0.5)
    plt.ylim(-0.5, msa.shape[0] - 0.5)
    plt.colorbar(label="Sequence identity to query", )
    plt.xlabel("Positions")
    plt.ylabel("Sequences")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}coverage_LDDT.png")
    ##################################################################
    plt.figure(figsize=(9, 6), dpi=300)
    ##################################################################
    plt.title(f"Predicted LDDT per position ({name})")
    for model_name, value in pae_plddt_per_model.items():
        plt.plot(value["plddt"], label=model_name)
    plt.ylim(0, 100)
    plt.ylabel("Predicted LDDT")
    plt.xlabel("Positions")
    plt.legend(loc="lower right")
    plt.savefig(f"{out_dir}/{name+('_' if name else '')}postion_LDDT.png")
    ##################################################################

    ##################################################################
    if is_multimer:
        for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
            # Create a new figure for each model
            fig = plt.figure(figsize=(9,6), dpi=300)
            # Set the title of the plot to the model name
            plt.title(f"Predicted allignment error ({name})")
            # Display the PAE as an image
            plt.imshow(value["pae"], label=model_name, cmap="Greens_r", vmin=0, vmax=30)
            # Add a colorbar to the plot
            plt.colorbar(label="Expected position error (Ångströms)")
            plt.ylabel("Aligned residue")
            plt.xlabel("Scored residue")
            # Save the plot as a PNG file, with the model name in the file name
            plt.savefig(f"{out_dir}/{name}_{model_name}_PAE.png")
            # Close the figure to free up memory
            plt.close(fig)
    ##################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir',dest='input_dir',required=True)
parser.add_argument('--name',dest='name')
parser.set_defaults(name='')
parser.add_argument('--output_dir',dest='output_dir')
parser.set_defaults(output_dir='')
args = parser.parse_args()

feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl','rb'))
is_multimer = ('result_model_1_multimer.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
# is_ptm = ('result_model_1_ptm.pkl' in [os.path.basename(f) for f in os.listdir(path=args.input_dir)])
# model_names = [f'{args.input_dir}/result_model_{f}{"_multimer" if is_multimer else "_ptm" if is_ptm else ""}.pkl' for f in range(1,6)]
model_names = sorted(glob.glob(f'{args.input_dir}/result_*.pkl'))

pae_plddt_per_model = get_pae_plddt(model_names, is_multimer=is_multimer)
generate_output_images(feature_dict, args.output_dir if args.output_dir else args.input_dir, args.name, pae_plddt_per_model, is_multimer)