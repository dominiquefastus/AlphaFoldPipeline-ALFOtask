#!/bin/bash
EDNA2_PATH=/data/staff/biomax/domfas/edna2_alphafold
. /home/domfas/Miniconda3/bin/activate
conda activate edna2       
export EDNA2_SITE=MAXIV_BIOMAX
export PATH=${EDNA2_PATH}/bin:$PATH
export PYTHONPATH=${EDNA2_PATH}/src