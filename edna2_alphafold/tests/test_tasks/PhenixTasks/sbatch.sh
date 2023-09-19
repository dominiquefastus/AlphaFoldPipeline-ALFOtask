#!/bin/bash
#SBATCH --job-name=Phenixtask
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=procpred_%j.out
#SBATCH --error=procpred_%j.err

#SBATCH --exclusive

source /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/PhenixTasks/start_sbatch.sh
cd  /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/PhenixTasks
python /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/PhenixTasks/ProcPredTask_exec_test.py
