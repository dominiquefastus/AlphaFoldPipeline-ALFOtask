#!/usr/bin/env bash
#SBATCH --job-name=AF_7QRZ
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

#SBATCH --exclusive
module purge
source /sw/tmp/z/m/activate
module add foss/2021b AlphaFold

export ALPHAFOLD_DATA_DIR=/sw/pkg/miv/mx/db/alphafold-2023a

export CWD=`pwd`
mkdir --parents alf_output/$SLURM_JOBID

cp /data/staff/biomax/domfas/data/7QRZ/7QRZ.fasta /local/slurmtmp.$SLURM_JOBID
cd /local/slurmtmp.$SLURM_JOBID

alphafold \
        --fasta_paths=7QRZ.fasta \
        --max_template_date=2022-01-01 \
        --db_preset=full_dbs \
        --model_preset=monomer \
        --output_dir=$CWD/alf_output/$SLURM_JOBID \
        --data_dir=$ALPHAFOLD_DATA_DIR

module purge
source /sw/tmp/z/m/deactivate