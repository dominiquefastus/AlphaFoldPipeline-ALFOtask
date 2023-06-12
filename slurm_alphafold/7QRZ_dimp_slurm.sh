#!/usr/bin/env bash
#SBATCH --job-name=AF_processing
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

module purge
module add gopresto CCP4/7.0.078-SHELX-ARP-8.0-17-PReSTO

cp 7qrz_phases.mtz alf_output/1408159/7QRZ
cd alf_output/1408159/7QRZ

dimple ranked_0_processed.pdb 7qrz_phases.mtz dimpleMR