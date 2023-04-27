#!/usr/bin/env bash
#SBATCH --job-name=AF_7QRZ
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

#SBATCH --exclusive
module purge
module add OpenSSL/1.0 fosscuda/2020b AlphaFold

export ALPHAFOLD_DATA_DIR=/sw/pkg/miv/mx/db/alphafold-2021b

export CWD=`pwd`
mkdir --parents alf_output/$SLURM_JOBID

cp /home/domfas/alphafold_pipeline/pipeline/7QRZ.fasta /local/slurmtmp.$SLURM_JOBID
cd /local/slurmtmp.$SLURM_JOBID

alphafold \
      --fasta_paths=7QRZ.fasta \
      --max_template_date=2020-05-14 \
      --output_dir=$CWD/alf_output/$SLURM_JOBID \
      --data_dir=$ALPHAFOLD_DATA_DIR