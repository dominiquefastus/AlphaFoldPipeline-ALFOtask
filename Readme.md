# AlphaFoldPipeline-ALFOtask (Implementation of a biomolecule structure solution pipeline based on AlphaFold prediction)

In this research project a standardised AlphaFold 2-based molecular replacement strategy was developed and implemented in an existing biomolecule structure solution pipeline at MAX IV Laboratory. The pipenline is designed to run on the MAX IV specific high perfomance cluster (HPC) which is build upon the LUNARC architecture. The biomolecular structure solution pipeline is programmed with the mind of user interaction at the beamline. Regardingly the pipeline is developed as a standalone program as well as implented in the on-site used edna2 framework. Therefore the project and repisotory is divided into two parts. The first part is the standalone program which is designed to run locally, but submit the individual computer extensive or program specific tasks to the MAX IV HPC. The second part is the implementation of the standalone program in the edna2 framework, which is directly run on the HPC and interacts with the MAX IV related computational infrastracture.

The standalone and edna2 implemented program are both designed to run the following tasks:
- AlphaFold 2 prediction (monomer or multimer)
- Processing of the predicted model
- Molecular replacement with the best model
- Refinement of the molecular replacement solution

The pipeline is also visualised in the following figure:

[![AlphaFold pipeline](figure)

## Requirements
To run the pipeline following resources and requirements are necessary to produce confident prediction:

### External programs installed on the HPC
- AlphaFold 2.1.1 
- Phenix 1.20.1
- CCP4 7.0.078

The programs are installed on the HPC and can be loaded with the module system and is therefore not necessary to be installed by the user. The module system is a software environment management system for the HPC. The commands to load the programs is explained in the sections.

### Dependency of gopresto or easybuild


## Standalone program

### Requirements
An access to the MAX IV HPC is needed, which is generally granted to users and scientists at MAX IV Laboratory. 

#### Python packages
- pathlib
- argparse
- subprocess
- logging
- shutil
- json
- time
- sys
- os

<br>
Self created utility functions:

- from utils import UtilsMonitor
- from utils import UtilsFileCheck


### 1. Run pipeline
The pipeline takes two input arguments, the fasta file of the protein of interest and the mtz file of the experimental data containing the information of the electron density map. 

Run the pipeline with the following command and example input files:
```
AlphaFold prediction pipeline

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA_PATH, --fasta FASTA_PATH
                        Path to fasta file to predict protein structure
  -m REFLECTIONDATA_PATH, --mtz REFLECTIONDATA_PATH
                        Mtz file to do molecular replacement
  -o OUTPUT, --output OUTPUT
                        Output directory for the results

Example: python3 AlphaFoldTask.py -f fasta_file -m mtz_file -o output_directory Author: D. Fastus
```

#### 1.1 Run the AlphaFold prediction
The first step creates a SLURM job and submits it to the MAX IV cluster. A fasta file for the protein of interest should be placed in the same directory as the 1_AF_prediction python script. The prediction is then automatically started for this protein fasta. The user will be reported when the job is done and the output file with all created predictions is completed.

The program will construct the following batch script and submit it to Slurm or the HPC:
```
#SBATCH --job-name=AF_7QRZ
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

#SBATCH --exclusive
module purge
module add fosscuda/2020b AlphaFold

export ALPHAFOLD_DATA_DIR=/db/alphafold-2021b

export CWD=`pwd`
mkdir --parents alf_output/$SLURM_JOBID

cp 7QRZ.fasta /local/slurmtmp.$SLURM_JOBID
cd /local/slurmtmp.$SLURM_JOBID

alphafold \
        --fasta_paths=7QRZ.fasta \
        --max_template_date=2022-01-01 \
        --db_preset=full_dbs \
        --model_preset=monomer \
        --output_dir=$CWD/alf_output/$SLURM_JOBID \
        --data_dir=$ALPHAFOLD_DATA_DIR
```

#### 1.2 Process the predicted model
From the previous AlphaFold 2 predicted structures, the model with the highest pLDDT (confidence for per-residue estimate) is used as a search model for molecular replacement. Before the model can be used for that further processing of the model needs to be done. The pLDDT values are replaced by pseudo b-values, as higher pLDDT values are associated with higher confidence in the prediction. This would ultimately lead to wrong b-value estimation, where higher values mean worse resolution. The processing will then trim low confident part of the model. Optionally the model is split into domains, which might enhance the molecular replacement solution for large proteins.

The program will construct the following batch script and submit it to Slurm or the HPC:
```
#SBATCH --job-name=AF_processing
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

module purge
module add gopresto Phenix/1.20.1-4487-Rosetta-3.10-PReSTO-8.0

cd alf_output/1478156/7QRZ
phenix.process_predicted_model ranked_0.pdb 
```

#### 1.3 Molecular replacement & refinement
For the molecular replacement and refinement the Dimpe pipeline from the CCP4 program is run automatically. 

The program will construct the following batch script and submit it to Slurm or the HPC:
```
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

dimple ranked_0_processed.pdb 7QRZ.mtz dimpleMR
```

The output of the pipeline is a refined model of the protein of interest and fitted into the electron density map. Since the AlphaFold yet produces different accurate models with rarely atomistic accurancy, further manual investigation of the model is needed. The model can be visualised with the CCP4 program Coot and refined here.

The ouptut folder named alf_output contains the following subfolders and files for the example monomer protein 7QRZ and the example multimer protein 8HUA:
```
├── 1478156
│   └── 7QRZ
│       ├── 1478161
│       ├── 1478162
│       ├── 7QRZ.mtz
│       ├── alphafold_1478156.err
│       ├── alphafold_1478156.out
│       ├── dimpleMR
│       │   ├── 01-rwcontents.log
│       │   ├── 02-phaser.log
│       │   ├── 03-unique.log
│       │   ├── 04-freerflag.log
│       │   ├── 05-cad.log
│       │   ├── 06-refmac5_jelly.log
│       │   ├── 07-refmac5_restr.log
│       │   ├── 08-find-blobs.log
│       │   ├── blob1-coot.py
│       │   ├── blob2-coot.py
│       │   ├── coot.sh
│       │   ├── dimple.log
│       │   ├── final.mmcif
│       │   ├── final.mtz
│       │   ├── final.pdb
│       │   ├── ini.pdb
│       │   ├── jelly.mmcif
│       │   ├── phaser.sol
│       │   ├── run-coot.py
│       │   ├── screen.log
│       │   └── workflow.pickle
│       ├── features.pkl
│       ├── msas
│       │   ├── bfd_uniclust_hits.a3m
│       │   ├── mgnify_hits.sto
│       │   ├── pdb_hits.hhr
│       │   └── uniref90_hits.sto
│       ├── ranked_0.pdb
│       ├── ranked_0_processed.pdb
│       ├── ranked_1.pdb
│       ├── ranked_2.pdb
│       ├── ranked_3.pdb
│       ├── ranked_4.pdb
│       ├── ranking_debug.json
│       ├── relaxed_model_1.pdb
│       ├── relaxed_model_2.pdb
│       ├── relaxed_model_3.pdb
│       ├── relaxed_model_4.pdb
│       ├── relaxed_model_5.pdb
│       ├── result_model_1.pkl
│       ├── result_model_2.pkl
│       ├── result_model_3.pkl
│       ├── result_model_4.pkl
│       ├── result_model_5.pkl
│       ├── timings.json
│       ├── unrelaxed_model_1.pdb
│       ├── unrelaxed_model_2.pdb
│       ├── unrelaxed_model_3.pdb
│       ├── unrelaxed_model_4.pdb
│       └── unrelaxed_model_5.pdb
└── 1478190
    └── 8HUA
        ├── 8HUA_raw.mtz
        ├── alphafold_1478190.err
        ├── alphafold_1478190.out
        ├── dimpleMR
        │   ├── 01-rwcontents.log
        │   ├── 02-truncate.log
        │   ├── 03-phaser.log
        │   ├── 04-unique.log
        │   ├── 05-freerflag.log
        │   ├── 06-cad.log
        │   ├── 07-refmac5_jelly.log
        │   ├── 08-refmac5_restr.log
        │   ├── 09-find-blobs.log
        │   ├── blob1-coot.py
        │   ├── blob2-coot.py
        │   ├── coot.sh
        │   ├── dimple.log
        │   ├── final.mmcif
        │   ├── final.mtz
        │   ├── final.pdb
        │   ├── ini.pdb
        │   ├── jelly.mmcif
        │   ├── phaser.sol
        │   ├── run-coot.py
        │   ├── screen.log
        │   └── workflow.pickle
        ├── features.pkl
        ├── msas
        │   ├── A
        │   │   ├── bfd_uniclust_hits.a3m
        │   │   ├── mgnify_hits.sto
        │   │   ├── pdb_hits.sto
        │   │   ├── uniprot_hits.sto
        │   │   └── uniref90_hits.sto
        │   ├── B
        │   │   ├── bfd_uniclust_hits.a3m
        │   │   ├── mgnify_hits.sto
        │   │   ├── pdb_hits.sto
        │   │   ├── uniprot_hits.sto
        │   │   └── uniref90_hits.sto
        │   ├── C
        │   │   ├── bfd_uniclust_hits.a3m
        │   │   ├── mgnify_hits.sto
        │   │   ├── pdb_hits.sto
        │   │   ├── uniprot_hits.sto
        │   │   └── uniref90_hits.sto
        │   └── chain_id_map.json
        ├── ranked_0.pdb
        ├── ranked_0_processed.pdb
        ├── ranked_1.pdb
        ├── ranked_1_processed.pdb
        ├── ranked_2.pdb
        ├── ranked_3.pdb
        ├── ranked_4.pdb
        ├── ranking_debug.json
        ├── relaxed_model_1_multimer.pdb
        ├── relaxed_model_2_multimer.pdb
        ├── relaxed_model_3_multimer.pdb
        ├── relaxed_model_4_multimer.pdb
        ├── relaxed_model_5_multimer.pdb
        ├── result_model_1_multimer.pkl
        ├── result_model_2_multimer.pkl
        ├── result_model_3_multimer.pkl
        ├── result_model_4_multimer.pkl
        ├── result_model_5_multimer.pkl
        ├── timings.json
        ├── unrelaxed_model_1_multimer.pdb
        ├── unrelaxed_model_2_multimer.pdb
        ├── unrelaxed_model_3_multimer.pdb
        ├── unrelaxed_model_4_multimer.pdb
        └── unrelaxed_model_5_multimer.pdb

```

## Edna2 implementation
The edna2 is a self developed framework for running multiple tasks related to macromolecular X-ray crystallography at synchrotron research facilities like the ESRF or MAX IV Laboratory.

edna2 needs to be seperately installed with using a conda environment.

### Requirements
An access to the MAX IV HPC is needed, which is generally granted to users and scientists at MAX IV Laboratory. The edna2 framework can also be run locally on a computer, but the program is designed to run on the HPC.

#### Python packages
- pathlib
- logging
- shutil
- json
- sys
- os

### Installation

### Execution of tasks
Despite the pipeline, the edna2 framework has defined tasks which are run independently and produce their own folder structure. An independent node needs to be set up to run the tasks on. This could be improved with more automation in the future like using the inbuild edna2 slurm system, which is not yet fully implemented. 

All input files for the specific task are defined in a json format as a design choice. All output information and parsing will be ultimately done with json files. The json files are also used to communicate between the different tasks. The json files are stored in the folder structure of the task. 

#### 2.1 Run the AlphaFold prediction task
The AlphaFold prediction task is identical to the standalone program. The task is run with the following commands and example input file:

JSON file for the example input file of the task to be executed:
```
{
    "fasta_path": "//data/7QRZ/7QRZ.fasta",
}
```

<br>

The slurm batch script is setup as following for the AlphaFold prediction task:
```
#!/bin/bash
#SBATCH --job-name=ALFOtask
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

#SBATCH --exclusive

source /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/AlphaFoldTask/start_sbatch.sh
cd  /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/AlphaFoldTask
python /data/staff/biomax/domfas/edna2_alphafold/tests/test_tasks/AlphaFoldTask/AlphaFold_exec_test.py
```

<br>

To submit or start the slurm batch script, the EDNA2 framework is setup with defining the path and activating the conda environment:
```
#!/bin/bash
EDNA2_PATH=/data/staff/biomax/domfas/edna2_alphafold
. /home/domfas/Miniconda3/bin/activate
conda activate edna2       
export EDNA2_SITE=MAXIV_BIOMAX
export PATH=${EDNA2_PATH}/bin:$PATH
export PYTHONPATH=${EDNA2_PATH}/src
```

#### 2.2 Run the processing task
The Process predicted model or procpred task is identical to the standalone program. The task is run within the Phenix task with the following commands and example input file:

JSON file for the example input file of the task to be executed:
```
{
    "PDB_file": "//data/alf_output/1468429/1JKB/ranked_0.pdb"
}
```

<br>

The slurm batch script is setup as following for the process predicted model task:
```
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

```

<br>
<br>

--- 

Please refer to the following sources:

[Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).](https://www.nature.com/articles/s41586-021-03819-2)

[Adams, Paul D et al. “PHENIX: a comprehensive Python-based system for macromolecular structure solution.” Acta crystallographica. Section D, Biological crystallography vol. 66,Pt 2 (2010): 213-21.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2815670/) 

