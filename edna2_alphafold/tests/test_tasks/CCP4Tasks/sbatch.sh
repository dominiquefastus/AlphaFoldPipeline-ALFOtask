#!/bin/bash
#SBATCH --job-name=Dimpletask
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=dimple_%j.out
#SBATCH --error=dimple_%j.err

#SBATCH --exclusive

source /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/CCP4Tasks/start_sbatch.sh
cd  /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/CCP4Tasks
python /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/CCP4Tasks/DimpleTask_exec_test.py
