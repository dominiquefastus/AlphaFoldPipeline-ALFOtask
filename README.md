# AlphaFoldPipeline-ALFOtask (Implementation of a biomolecule structure solution pipeline based on AlphaFold prediction)

```
Author:     Dominique Fastus
Date:       2021-10-01

Github:     https://github.com/dominiquefastus/AlphaFoldPipeline-ALFOtask
Tested:     MAX IV Cluster, LUNARC architecture (stable)
```
In this research project a standardised AlphaFold 2-based molecular replacement strategy was developed and implemented in an existing biomolecule structure solution pipeline at MAX IV Laboratory. The pipeline is designed to run on the MAX IV specific high perfomance cluster (HPC) which is build upon the LUNARC architecture. The biomolecular structure solution pipeline is programmed with the mind of user interaction. Regardingly, the pipeline is developed as a standalone program, as well as implented in the MAXIV used EDNA2 framework. Therefore, the project and repository is divided into two parts. The first part is the standalone program, which is designed to run locally, but submit the individual computer extensive or program specific tasks to the MAX IV HPC. The second part is the implementation of the standalone program in the EDNA2 framework, which is directly run on the HPC and interacts with the MAX IV and beamline related computational infrastracture.

The scripts also provide general plotting for molecular replacement results of different metrics and tools to interact and manipulate pdb files that are useful outside of the pipeline.

The standalone and EDNA2 implemented pipelines are both designed to run the following tasks:
- AlphaFold 2 prediction (monomer or multimer)
- Processing of the predicted model
- Molecular replacement with the best model
- Refinement of the molecular replacement solution


## Requirements
To run the pipeline, following resources and requirements are necessary to produce the solutions:

### External programs installed on the HPC

MAX IV Computer Cluster:
```
- AlphaFold             2.1.1 
- Phenix                1.20.1
- CCP4                  7.0.078
```

The programs are installed on the HPC and can be loaded with the module system and are therefore not necessary to be installed by the user. The module system is a software environment management system for the HPC. The commands to load the programs are explained in the sections.  

### Dependency on PReSTO and EasyBuild
"PReSTO is a software stack for integrated structural biology adapted to high performance computing resources at the National Academic Infrastructure for Su-
percomputing in Sweden (NAISS) and the local MAX IV compute cluster" [1]. PReSTO is used in the project to load the modules Phenix and CCP4. The exact loading of the software with PReSTO is explained in the sections.

Easybuild [2] is a frameworl to build and install software on the HPC. Easybuild is used in the project to load the module AlphaFold. The exact loading of the software with easybuild is explained in the sections. The easybuild framework implemented AlphaFold version 2.1.1, seemed to be the most stable and reliable build of AlphaFold on the cluster, but this should updated to the PReSTO version in the future.

1 PReSTO_docs_2023_Tallberg.pdf. (n.d.). Retrieved November 1, 2023, from https://www.nsc.liu.se/support/presto/PReSTO_docs_2023_Tallberg.pdf
<br>2 EasyBuild—Building software with ease. (n.d.). Retrieved November 1, 2023, from https://docs.easybuild.io/


### Slurm workload manager
Slurm is a workload cluster management and job scheduling system, which is used in the project to submit and monitor the individual computing tasks to the HPC. This eases the parallel run of multiple jobs or tasks and allows to better destribute the computing ressources to multiple users.

## Standalone program

### Requirements
An access to the MAX IV HPC is needed, which is generally granted to users and scientists at MAX IV Laboratory. 

#### Python packages
The environment uses Python =3.8.15. Only standard python packages and self created utility functions are used in the standalone program:
```python
- from utils import UtilsMonitor
- from utils import UtilsFileCheck
```

### 1. Run pipeline
The pipeline takes two input arguments, the fasta file of the protein of interest and the mtz file of the experimental data containing the information of the electron density map. 

Run the pipeline with the following command and example input files:
```
usage: AlphaFoldTask.py [-h] -f FASTA_PATH -m REFLECTIONDATA_PATH [-o OUTPUT] [-p PDB_PATH] [-jp] [-pp] [-t]

AlphaFold prediction pipeline

options:
  -h, --help            show this help message and exit
  -f FASTA_PATH, --fasta FASTA_PATH
                        Path to fasta file to predict protein structure
  -m REFLECTIONDATA_PATH, --mtz REFLECTIONDATA_PATH
                        Mtz file to do molecular replacement
  -o OUTPUT, --output OUTPUT
                        Output directory for the results
  -p PDB_PATH, --pdb PDB_PATH
                        Path to pdb file to predict protein structure
  -jp, --just-predict   Only predict the structure, no molecular replacement
  -pp, --predict-process
                        Predict the structure and process it, no molecular replacement
  -t, --tidy            Tidy up temporary files

Example: python3 AlphaFoldTask.py -f fasta_file -m mtz_file -o output_directory
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

The AlphaFold job script will be automatically adapted based on the number of sequences in the fasta input file, which indicates monomer or multimer proteins. So only `--model_preset=multimer` will be changed accordingly.

#### 1.2 Process the predicted model
From the previous AlphaFold 2 predicted structures, the model with the highest pLDDT (confidence for per-residue estimate) is used as a search model for molecular replacement. Before the model can be used for that, further processing of the model needs to be done. The processing will also trim predicted parts with low confidence (pLDDT <70) of the model. 

The processing is based on the method suggested by Oeffner et all [3] and was implemented with python into the pipeline. Whereas proviously realying on the Phenix program, the processing is now done with python calculations. The implemented processing lacks the possibility to split the model into domains, thus the phenix program is still provided if this is needed.

To run the processing with Phenix, the pipeline would construct the following batch script and submit it to Slurm or the HPC:
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

3. Oeffner RD, Croll TI, Millán C, Poon BK, Schlicksup CJ, Read RJ, et al. Putting AlphaFold models to work with phenix.process_predicted_model and ISOLDE. Acta Cryst D. 2022 Nov 1;78(11):1303–14. 

For the basic processing of the inbuild version of phenix.process_predicted_model, the implemented version without domain splitting is enough to archieve good search models for phasing and is prefered.

#### 1.3 Molecular replacement & refinement
For the molecular replacement and refinement the Dimple pipeline from the CCP4 program is run automatically. 

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

cp /data/staff/biomax/user/data/7QRZ/7QRZ.mtz alf_output/1478156/7QRZ
cd alf_output/1478156/7QRZ

dimple ranked_0_processed.pdb 7QRZ.mtz dimpleMR
```
<br>
The output of the pipeline is a refined model of the protein of interest based on the fitting to the electron density map. Since the AlphaFold produces different accurate models with varying atomistic accurancy, further manual investigation of parts or the full model is needed. 
<br>
<br>
The directory is structured with the job id on the first level and with another subfolder indicating the protein name. It contains the input file for the pipeline run, the logging error files as well as the predicted models in pdb format and additional files from the alphafold prediction, which are explained in more detail at https://github.com/google-deepmind/alphafold. The directory also includes the dimpleMR folder with log files, the final refined model in pdb format and the corresponding reflection data in mtz format after molecular replacement and refinement.

## EDNA2 implementation
EDNA2 [4] is a self developed framework for running multiple tasks related to macromolecular X-ray crystallography at synchrotron research facilities like the ESRF or MAX IV Laboratory.

EDNA2 needs to be seperately installed with using a conda environment. The environment uses Python =3.8.15.

[4] Incardona, M.-F., Bourenkov, G. P., Levik, K., Pieritz, R. A., Popov, A. N., & Svensson, O. (2009). EDNA: A framework for plugin-based applications applied to X-ray experiment online data analysis. Journal of Synchrotron Radiation, 16(6), 872–879. https://doi.org/10.1107/S0909049509036681

### Requirements
An access to the MAX IV HPC is needed, which is generally granted to users and scientists at MAX IV Laboratory. The EDNA2 framework can also be run locally on a computer, but the program is designed to run on the HPC. 

#### Python packages
EDNA2 requires several python packages to be installed in the conda environment. The packages are listed in the edna2-alphafold/requirements.txt file or below:
```
- cctbx-base                2023.5           py38h81fdc0f_0    conda-forge
- graypy                    2.1.0                    pypi_0    pypi
- h5py                      3.7.0                    pypi_0    pypi
- hdf5plugin                4.1.1                    pypi_0    pypi
- matplotlib                3.6.3                    pypi_0    pypi
- xmltodict                 0.13.0                   pypi_0    pypi
- suds-jurko                0.6           py38h578d9bd_1005    conda-forge
- jsonschema                4.17.3                   pypi_0    pypi
- numpy                     1.24.3           py38hf6e8229_1 
- fabio                     2022.12.1                pypi_0    pypi
- requests                  2.28.2                   pypi_0    pypi
- distro                    1.8.0                    pypi_0    pypi
- scipy                     1.10.0                   pypi_0    pypi
- billiard                  4.1.0                    pypi_0    pypi
```


### Installation
To install the EDNA2 framework, a conda environment need to be setup:

- Download and install miniconda:
```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-py38_23.3.1-0-Linux-x86_64.sh
$ bash Miniconda3-latest-Linux-x86_64.sh
```

- Create a conda environment with specific python version:
```
(base) $ conda create -n EDNA2 python=3.8
```

- Activate the conda environment, noticable by the name in front of the command line:
```
(base) $ conda activate EDNA2
(EDNA2) $
```

- Clone the EDNA2 repository from the github:
```
(EDNA2) $ git clone git@github.com:aaronfinke/EDNA2.git
(EDNA2) $ cd EDNA2
```

Use the EDNA2_alphafold branch which contaions the AlphaFold implementation and a yaml setup file.
- Install necessary packages:
```
(EDNA2) $ conda install -c conda-forge --file requirements.txt
```

- After that EDNA2 can be installed with the setup.py file in this environment:
```
(EDNA2) $ cd /path/to/EDNA2
(EDNA2) $ python setup.py install
```
- Alternatively, the EDNA2 environment can be created with the environment.yml file:
```
(base) $ conda env create -f environment.yml
(EDNA2) $ cd /path/to/EDNA2
(EDNA2) $ python setup.py install
```

### Execution of tasks
Despite the pipeline, the EDNA2 framework has defined tasks which are run independently and produce their own folder structure. An independent node needs to be set up to run the tasks on. This could be improved with more automation in the future like using the inbuild EDNA2 slurm system, which is not yet fully implemented. 

All input files for the specific task are defined in a json format as a design choice. All output information and parsing will be ultimately done with json files. The json files are also used to communicate between the different tasks. The json files are stored in the folder structure of the task. 

<br>

To submit or start the specific task, the EDNA2 framework is setup each time with defining the path and activating the conda environment [file name: start_sbatch.sh]:
```
#!/bin/bash
EDNA2_PATH=/data/staff/biomax/user/EDNA2_alphafold
. /home/user/Miniconda3/bin/activate
conda activate EDNA2       
export EDNA2_SITE=MAXIV_BIOMAX
export PATH=${EDNA2_PATH}/bin:$PATH
export PYTHONPATH=${EDNA2_PATH}/src
```


#### 2.1 Run the AlphaFold prediction task
The AlphaFold prediction task is identical to the standalone program. The task is run with the following commands and example input file:

JSON file for the example input file of the task to be executed:
```
{
    "fasta_path": "//data/7QRZ/7QRZ.fasta",
}
```

<br>

The slurm batch script is setup as following for the AlphaFold prediction task [file name: alf_sbatch.sh] by sourcing the start_sbatch.sh file and running the python test execution script:
```
#!/bin/bash
#SBATCH --job-name=ALFOtask
#SBATCH --partition=v100
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=alphafold_%j.out
#SBATCH --error=alphafold_%j.err

#SBATCH --exclusive

source /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/AlphaFoldTask/start_sbatch.sh
cd  /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/AlphaFoldTask
python /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/AlphaFoldTask/AlphaFold_exec_test.py
```

#### 2.2 Run the processing task
The Process predicted model or procpred task is run with the phenix program and does not provide the faster, simpler inbuild python version yet. The task is run within the Phenix task with the following commands and example input file:

JSON file for the example input file of the task to be executed:
```
{
    "PDB_file": "//data/alf_output/1468429/1JKB/ranked_0.pdb"
}
```

<br>

The slurm batch script is setup as following for the process predicted model task [file name: ppm_sbatch.sh] by sourcing the start_sbatch.sh file and running the python test execution script:
```
#!/bin/bash
#SBATCH --job-name=Phenixtask
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=procpred_%j.out
#SBATCH --error=procpred_%j.err

#SBATCH --exclusive

source /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/PhenixTasks/start_sbatch.sh
cd  /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/PhenixTasks
python /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/PhenixTasks/ProcPredTask_exec_test.py

```

#### 2.2 Run the molecular replacement task
The Dimple task is run with the CCP4 program and does not provide the faster, simpler inbuild python version yet. The task is run within the CCP4 task with the following commands and example input file:

JSON file for the example input file of the task to be executed:
```
{
    "PDB_file": "//data/alf_output/1468429/1JKB/ranked_0_processed.pdb",
    "MTZ_file": "/data//data/7QRZ/7QRZ.mtz"
}
```

The slurm batch script is setup as following for the dimple task [file name: dim_sbatch.sh] by sourcing the start_sbatch.sh file and running the python test execution script:
```
#!/bin/bash
#SBATCH --job-name=Dimpletask
#SBATCH --mem=0
#SBATCH --time=01-00:00
#SBATCH --output=dimple_%j.out
#SBATCH --error=dimple_%j.err

#SBATCH --exclusive

source /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/CCP4Tasks/start_sbatch.sh
cd  /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/CCP4Tasks
python /data/staff/biomax/user/EDNA2_alphafold/tests/test_tasks/CCP4Tasks/DimpleTask_exec_test.py
```

The actual test tasks is then run like following:
```
sbatch [alf_/ppm_/dim_]sbatch.sh
```
## Utilsplots for visualisation
### Requirements
Results from the standalone or EDNA2 implemented pipeline.

#### Python packages
The environment uses Python >=3.8. The following packages are used in the plotting scripts:
```
- matplotlib                3.6.3                    pypi_0    pypi
- biopandas                 0.4.1                    pypi_0    pypi
- biopython                 1.81             py38h1de0b5d_0    conda-forge
- pandas                    2.0.3            py38h1128e8f_0 
- numpy                     1.24.3           py38hf6e8229_1  
```
#### Run the plotting
- To plot the b-factors along the residue positions of the protein, the following command can be used:
```
usage: bfactor_plot.py [-h] [-l LABELS [LABELS ...]] [-a] pdb_files [pdb_files ...]

Plot B-factors for multiple PDB files.

positional arguments:
  pdb_files             Paths to PDB files

optional arguments:
  -h, --help            show this help message and exit
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...]
                        Labels for the PDB files
  -a, --align           Align residues across different PDB files
```
- To plot the boxplot of the plDDT values of the AlphaFold prediction against the b-factors of the reference models, the following command can be used:
```
usage: plddt_bval.py [-h] -p ALPHAFOLD_PDB REFERENCE_PDB [-l LABEL]

Process AlphaFold and reference PDB files.

optional arguments:
  -h, --help            show this help message and exit
  -p ALPHAFOLD_PDB REFERENCE_PDB, --pair ALPHAFOLD_PDB REFERENCE_PDB
                        Specify a pair of PDB files
  -l LABEL, --label LABEL
                        Label for each pair
```
- To plot the rfree, LLG and TFZ values of the molecular replacement results, the following command can be used:
```
usage: MR_metrics.py [-h] logfile

Plot data from dimple.log file.

positional arguments:
  logfile     Path to the dimple.log file
```
- To plot the overall metric statistics of the prediction and molecular replacement results, the following command can be used:
```
usage: overall_MR.py [-h] csv_file

Generate protein metrics plots from CSV data.

positional arguments:
  csv_file    Path to the CSV file containing protein metrics.

optional arguments:
  -h, --help  show this help message and exit
```
#### Supplementary scripts
- To plot the AlphaFold pLDDT for the residues and the PAE for multimers an adaptet script from AlphaFold and friends on the HPC. (n.d.). Retrieved May 6, 2024, from https://elearning.vib.be/courses/alphafold/ (https://github.com/jasperzuallaert/VIBFold/blob/main/visualize_alphafold_results.py) is used with the following command:
```
usage: prediction_plot.py [-h] --input_dir INPUT_DIR [--name NAME] [--output_dir OUTPUT_DIR]

options:
  -h, --help            show this help message and exit
  --input_dir INPUT_DIR
  --name NAME
  --output_dir OUTPUT_DIR
  ```
<br>
<br>

--- 

Please refer to the following sources:

[Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).](https://www.nature.com/articles/s41586-021-03819-2)

[Adams, Paul D et al. “PHENIX: a comprehensive Python-based system for macromolecular structure solution.” Acta crystallographica. Section D, Biological crystallography vol. 66,Pt 2 (2010): 213-21.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2815670/) 

