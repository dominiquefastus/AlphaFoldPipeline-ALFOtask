#!/usr/bin/env python3

"""
Script for plotting the bvalue per residue

Author:     D. Fastus
"""
#%%
from biopandas.pdb import PandasPdb

import matplotlib.pyplot as plt
from matplotlib import style


# fetch reference pdb from rcsb
reference_pdb = PandasPdb().fetch_pdb("7qrz")

# or fetch structure from alphafold database
# PandasPdb().fetch_pdb(uniprot_id='Q5VSL9', source="alphafold2-v2")

# load in the pdb to analyze and compare with like (predicted_model, processed_model, refined_model)
predicted_model = PandasPdb().read_pdb(path="/data/staff/biomax/domfas/AlphaFold_project/slurm_alphafold/alf_output/1478156/7QRZ/ranked_0.pdb")
processed_model = PandasPdb().read_pdb(path="/data/staff/biomax/domfas/AlphaFold_project/slurm_alphafold/alf_output/1478156/7QRZ/ranked_0_processed.pdb")
refined_model = PandasPdb().read_pdb(path="/data/staff/biomax/domfas/AlphaFold_project/slurm_alphafold/alf_output/1478156/7QRZ/dimpleMR/final.pdb")

refined_model.df['ATOM']['b_factor'].plot(kind='line', label='predicted model')
reference_pdb.df['ATOM']['b_factor'].plot(kind='line', label='deposited model')
plt.xlim((0,1400))
plt.legend(loc='upper right')
style.use('ggplot')
plt.title('Distribution of B-Factors')
plt.xlabel('B-factor')
plt.ylabel('count')
plt.show()


# %%