""" 
This is the script to produce AlphaFold prediction quality metrics plots.

The general code for the plots was adapted and modified from:
https://github.com/jasperzuallaert/VIBFold/blob/main/visualize_alphafold_results.py

The adapted code is also featured in:
https://blog.biostrand.be/explained-how-to-plot-the-prediction-quality-metrics-with-alphafold2

"""

#%%
import os
import glob
import math
import json
import pickle
import matplotlib.pyplot as plt
import numpy as np

folder = "/data/staff/biomax/domfas/AlphaFold_project/slurm_alphafold/alf_output/1530559/5NJM"
name = "5NJM AlphaFold prediction"

class ARG:
    def __init__(self, folder):
        self.input_dir = folder
        self.output_dir = folder
        self.name = folder

args = ARG(folder)

# open the rankings of the predicted models as json files
with open(os.path.join(folder, "ranking_debug.json"), 'r') as f:
    ranking_dict = json.load(f)

feature_dict = pickle.load(open(f'{args.input_dir}/features.pkl','rb'))

if 'result_model_1_multimer.pkl' in [os.path.basename(f) for f in os.listdir(path=folder)]:
    is_multimer = True
else:
    is_multimer = False
    

if is_multimer == False:
    model_dicts = [pickle.load(open(f'{args.input_dir}/result_model_{f}{"_multimer" if is_multimer else ""}{"_ptm" if is_multimer else ""}.pkl','rb')) for f in range(1,6)]
else:
    model_dicts = [pickle.load(open(f'{args.input_dir}/result_model_{f}{"_multimer" if is_multimer else ""}{"_ptm" if is_multimer==False else ""}.pkl','rb')) for f in range(1,6) for g in range(5)]

    
# define function to get the plddt score
def get_pae_plddt(model_dicts):
    out = {}
    
    for i,d in enumerate(model_dicts):
        out[f'model_{i+1}'] = {'plddt': d['plddt'],
                               'pae':d['distogram']}
    return out

def generate_MSA_image(feature_dict, name):
            msa = feature_dict['msa']
            seqid = (np.array(msa[0] == msa).mean(-1))
            seqid_sort = seqid.argsort()
            non_gaps = (msa != 21).astype(float)
            non_gaps[non_gaps == 0] = np.nan
            final = non_gaps[seqid_sort] * seqid[seqid_sort, None]
            
            # plot MSA with coverage
            plt.figure(figsize=(10, 6), dpi=100)
            plt.title(f"Sequence coverage for ({name})")
            plt.imshow(final, interpolation='nearest', aspect='auto',
                       cmap="rainbow_r", vmin=0, vmax=1, origin='lower')
            plt.plot((msa != 21).sum(0), color='black')
            plt.xlim(-0.5, msa.shape[1] - 0.5)
            plt.ylim(-0.5, msa.shape[0] - 0.5)
            plt.colorbar(label="Sequence identity to query", )
            plt.xlabel("Positions")
            plt.ylabel("Sequences")
            
def generate_plddt_image(feature_dict, model_dicts, ranking_dict, name):
            msa = feature_dict['msa']
            seqid = (np.array(msa[0] == msa).mean(-1))
            seqid_sort = seqid.argsort()
            non_gaps = (msa != 21).astype(float)
            non_gaps[non_gaps == 0] = np.nan
            final = non_gaps[seqid_sort] * seqid[seqid_sort, None]
            
            pae_plddt_per_model = get_pae_plddt(model_dicts)
            
            # plot plddt per position
            plt.figure(figsize=(10, 6), dpi=100)
            plt.title(f"Predicted LDDT per position ({name})")
            
            s = 0
            for model_name, value in pae_plddt_per_model.items():
                if is_multimer:
                    plt.plot(value["plddt"], label=f"{model_name}_multimer, plddts: {round(list(ranking_dict['iptm+ptm'].values())[s], 6)}")
                else:
                   plt.plot(value["plddt"], label=f"{model_name}, plddts: {round(list(ranking_dict['plddts'].values())[s], 6)}")
                    
                s += 1
                plt.legend()
                plt.ylim(0, 100)
                plt.ylabel("Predicted LDDT")
                plt.xlabel("Positions")

'''        
def generate_pae_image(model_dicts):
    model_names = sorted(glob.glob(f'{folder}/result_*.pkl'))
    pae_plddt_per_model = get_pae_plddt(model_dicts)
    
    num_models = 5 # columns
    num_runs_per_model = math.ceil(len(model_names)/num_models)
    fig = plt.figure(figsize=(3 * num_models, 2 * num_runs_per_model), dpi=100)
    for n, (model_name, value) in enumerate(pae_plddt_per_model.items()):
        print(value)
        plt.subplot(num_runs_per_model, num_models, n + 1)
        plt.title(model_name)
        plt.imshow(value["pae"], label=model_name, cmap="bwr", vmin=0, vmax=30)
        plt.colorbar()
    fig.tight_layout()
'''      

#%%
generate_MSA_image(feature_dict, name)




# %%
generate_plddt_image(feature_dict, model_dicts, ranking_dict, name)




# %%
