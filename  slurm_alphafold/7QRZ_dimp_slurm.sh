#!/usr/bin/env bash
#SBATCH --job-name=AF_processing
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

module purge
module add gopresto CCP4/7.0.078-SHELX-ARP-8.0-17-PReSTO

cp /data/staff/biomax/domfas/data/7QRZ/7QRZ.mtz alf_output/1478156/7QRZ
cd alf_output/1478156/7QRZ

dimple ranked_0_processed.pdb /data/staff/biomax/domfas/data/7QRZ/7QRZ.mtz dimpleMR