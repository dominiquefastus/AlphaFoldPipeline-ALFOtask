#!/usr/bin/env bash
#SBATCH --job-name=AF_processing
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

module purge
module add gopresto Phenix/1.20.1-4487-Rosetta-3.10-PReSTO-8.0

cd alf_output/1464798/7QRZ
phenix.process_predicted_model ranked_0.pdb 