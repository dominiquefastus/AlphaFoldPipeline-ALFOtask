#!/usr/bin/env python3

"""
This is a automated python based SLURM job submitting pipeline for AlphaFold prediction. 

It serves as a wrapper for the AlphaFold prediction and uses the installed AlphaFold module on the server.
It allows to predict the structure of a protein from a fasta file and do molecular replacement and refinement with DIMPLE.

Users can provide a fasta file and a mtz file to do molecular replacement and refinement with DIMPLE. 
All classes are also implemented as seperat tasks in the edna framework and can be used as such.

The pipeline is divided in three steps represented by three classes:
- class ALFOpred: initiate the alphafold prediction as a slurm job by using the installed alphafold module on the server
- class procALFO: process the best model from the prediction with phenix.process_predicted_model or implemented version of it
- class mrALFO: do molecular replacement and refinement with dimple

To run the pipeline following arguments are necessary or optional:
usage: AlphaFoldTask.py [-h] -f FASTA_PATH -m REFLECTIONDATA_PATH [-o OUTPUT] [-p PDB_PATH] [-jp] [-pp] [-t]

AlphaFold prediction pipeline

optional arguments:
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

Author:     D. Fastus
"""

# necessary imports for the script
from pathlib import Path
import argparse
import logging
import shutil
import math
import json
import sys
import re
import os

# import module to monitor the SLURM jobs
from utils import UtilsMonitor

# import module to check the output files
from utils import UtilsFileCheck

# build a logger for the SLURM script
logging.basicConfig(filename=f'job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# add a handler to print the log messages to stdout
# also define the format of the log messages
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


# Add argparse to the script to get diffent input files
parser = argparse.ArgumentParser(
                    prog="AlphaFoldTask.py",
                    description="AlphaFold prediction pipeline",
                    epilog='''Example: python3 AlphaFoldTask.py -f fasta_file -m mtz_file -o output_directory\n\n
                    Author: D. Fastus''')

# Add arguments
parser.add_argument("-f", "--fasta", dest="fasta_path", required=True, 
                    help="Path to fasta file to predict protein structure")
parser.add_argument("-m", "--mtz", dest="reflectionData_path", required=True,
                    help="Mtz file to do molecular replacement")
parser.add_argument("-o", "--output", default="alf_output", required=False, 
                    help="Output directory for the results")

# alternatve provide pdb and mtz file
parser.add_argument("-p", "--pdb", dest="pdb_path", required=False, 
                    help="Path to pdb file to predict protein structure")

# option to just predict or predict and process the model
parser.add_argument("-jp", "--just-predict", action="store_true", dest="just_predict", required=False, 
                   help="Only predict the structure, no molecular replacement")
parser.add_argument("-pp", "--predict-process", action="store_true", dest="predict_process", required=False, 
                   help="Predict the structure and process it, no molecular replacement")

# option to tidy up temporary files
parser.add_argument("-t", "--tidy", action="store_true", dest="tidy", required=False, 
                   help="Tidy up temporary files")

# Parse arguments  
args = parser.parse_args()

# the pipeline is designed in a object oriented way to keep the code clean and readable but als keep the structure 
# of the edna framework for easier implementation 

# alphafold is a deep learning based method for protein structure prediction and available through: 
# https://github.com/google-deepmind/alphafold 
class ALFOpred():
    """
    Initiates a AlphaFold prediction SLURM job and automatically sets the monomer or multimer preset
    """

    # base id where everything is stored
    def __init__(self):
        self.job_id = None

    # create sbatch file for server to run alphafold with monomer or multimer preset
    # uses the most stable and reliable version of alphafold on the server: AlphaFold/2.1.1
    # the methods of the class takes parameters to create the sbatch file and initiate the job 
    # to set up the node at the cluster specific requirements are giving
    def sbatch_AFpred(self, preset, jobName, fastaName, fastaPath, partition = "v100", mem = 0, time = "01-00:00"):
        # create the script for the sbatch file
        script = "#!/usr/bin/env bash\n"

        # setting up the cluster environment or job specifications
        script += f"#SBATCH --job-name=AF_{fastaName}\n"

        # partition is optional, if not given the default partition will be used
        if partition is not None:
            script += f"#SBATCH --partition={partition}\n"
        
        # define the memory and time for the job at the node of the cluster
        script += f"#SBATCH --mem={mem}\n"
        script += f"#SBATCH --time={time}\n"

        # create a output and error file for the job
        script += "#SBATCH --output=alphafold_%j.out\n"
        script += "#SBATCH --error=alphafold_%j.err\n\n"

        # make the node exclusive for the job for better performance
        script += "#SBATCH --exclusive\n"

        # alphafold commands and requirements
        # load the alphafold module and purge the current modules
        script += "module purge\n"   
        script += "module add fosscuda/2020b AlphaFold\n\n"

        # define alphafold database directory for reference to run the msas
        script += "export ALPHAFOLD_DATA_DIR=/sw/pkg/miv/mx/db/alphafold-2021b\n\n"

        # define the working directory and create a output directory for the job
        script += "export CWD=`pwd`\n"
        script += f"mkdir --parents alf_output/$SLURM_JOBID\n\n"

        # copy the fasta file to the output directory and change the directory to the output directory where the prediction will be done
        script += f"cp {fastaPath} /local/slurmtmp.$SLURM_JOBID\n"
        script += f"cd /local/slurmtmp.$SLURM_JOBID\n\n"

        # run alphafold
        # change accordingly the model_preset: monomer, monomer_casp14 or monomer_ptm 
        # change accordingly the db_preset: full_dbs or reduced_dbs
        # for faster prediction use the reduced_dbs and monomer_ptm
        basename = os.path.basename(fastaPath)
        script += f"""alphafold \\
        --fasta_paths={basename} \\
        --max_template_date=2022-01-01 \\
        --db_preset=full_dbs \\
        --model_preset={preset} \\
        --output_dir=$CWD/alf_output/$SLURM_JOBID \\
        --data_dir=$ALPHAFOLD_DATA_DIR"""
        
        # name the sbatch file and write the script to the file
        shellFile = f"{jobName}_slurm.sh"

        # wrtie the script to the file 
        with open(shellFile,"w") as file:
            file.write(script)

    # methods to check if the output are complete and have the right format
    # only the most important files are checked
    def check_out(self, out_dir, preset):
        # List of files to check for in monomer prediction
        files_to_check_mono = ['ranked_0.pdb', 'relaxed_model_1.pdb', 'result_model_1.pkl', 'unrelaxed_model_1.pdb', "ranking_debug.json"]

        # List of files to check for in multimer prediction
        files_to_check_multi = ['ranked_0.pdb', 'relaxed_model_1_multimer.pdb', 'result_model_1_multimer.pkl', 'unrelaxed_model_1_multimer.pdb', "ranking_debug.json"]

        # Check if all defined files are present in the directory
        # use the method loop_files from the UtilsFileCheck module
        # change for monomer or multimer prediction
        if preset == "monomer":
            UtilsFileCheck.loop_files(output_dir=out_dir, files_to_check=files_to_check_mono)
        else:
            UtilsFileCheck.loop_files(output_dir=out_dir, files_to_check=files_to_check_multi)

        # if all files are present the script will continue and log a message  
        logger.info(f"The files are successfully generated in {out_dir}")
        return True

    # method to run the prediction as a slurm job
    # takes the fasta file as an argument (args), the only mandatory argument for the aalphafold prediction
    def run(self, args):
        # check if it's a fasta file and get the name for the job
        try:
            with open(args,mode="r") as file:
                line = file.read()

                # reads the first line of the fasta file and check if it starts with a '>'
                # if not the script will stop and give an error message, which will be logged
                if not line.startswith('>'):
                    logger.error("The input is not a fasta file!")
                    sys.exit(1)

                # if the first line starts with a '>', it will be counted how many '>' or sequences are in the file
                # by that the script can determine if it's a monomer or multimer prediction and c hange the preset accordingly
                else:
                    if line.count('>') == 1:
                        preset = "monomer"
                    else:
                        preset = "multimer"
                    # get the name of the fasta file without the extension
                    line = line.strip()
                    fasta_name = line[1:5]

        # if there is an error with the fasta file the script will stop and give an error message, which will be logged
        # most likely the file doesn't exist or can't be opened
        except Exception as e:
            logger.error(f"{args} can not be open or does not exist!", exc_info=True)
            sys.exit(1)
        
        # create the sbatch file depending on the fasta file
        # the preset is determined by the number of sequences in the fasta file
        # get from the fasta file the name of the protein
        # instantiate the class and run the sbatch file an run the prediction
        ALFOpred.sbatch_AFpred(preset=preset, jobName=f"AFpred_{fasta_name}", fastaName=fasta_name, fastaPath=args)
    
        # get the job id from the monitor and print the status of the job to stdout
        # this can be replaced by a email notification and using --wait in the sbatch file
        job_id = UtilsMonitor.monitor_job(script=f"AFpred_{fasta_name}_slurm.sh", name=f"Predicting ({fasta_name})")
        self.job_id = job_id

        # get the output directory to check if the prediction was successful and move some files
        output_dir = f"alf_output/{job_id}/{fasta_name}"

        # check if the prediction was successful
        # check if alphafold output is complete (files are present)
        ALFOpred.check_out(out_dir=output_dir, preset=preset)

        # move the slurm output files to the output directory
        move_file = ["job.log", f"alphafold_{job_id}.err", f"alphafold_{job_id}.out"]
        for file in move_file:
            shutil.move(file, output_dir)

        # if the prediction was successful the script will continue and log a message
        logger.info("The prediction was successful!\n\n")


        # return the output directory and the fasta name for the next step
        return output_dir, fasta_name, job_id, preset


# the model from the prediction is processed and the best model is selected
# for further molecular replacement and refinement the model is processed with phenix
# phenix is a software for the automated determination of molecular structures using X-ray crystallography and other methods
# the method process_predicted_model is used to process the model and is available through: https://www.phenix-online.org/documentation/reference/process_predicted_model.html
# basically the method replaces the B-factor values with estimated B values and removes low-confidence residues and (optional splits into domains, not used here)
class procALFO():
    """
    AlphaFold predicted model selection and processing with phenix.process_predicted_model or implemented version of it
    """

    def __init__(self, alphafold_job_id):
        self.alphafold_job_id = alphafold_job_id

    # method to choose the best model from the prediction
    # the best model is chosen by the plddts score
    # this is only necessary for loging purposes and to give the user information about the best model
    # alphafold ranks the models automatically from best to worst and the best model is always ranked as 0
    def get_model(self, infile):
        try:
            """
            # open the ranking_debug.json file and load the json file
            json_file = os.path.join(infile, "ranking_debug.json")
            file = open(json_file)
            ranking = json.load(file)

            # give information about the best model
            best_model = ranking["order"][0]

            # monomer and multimer give different scoring units for the best model
            # need to try both units otherwise the script will fail
            try:
                plddts = ranking["plddts"][best_model]
                logger.info(f"The best predicted model by AlphaFold is {best_model} with a plddts of {plddts}")
            except KeyError:
                iptm_ptm = ranking["iptm+ptm"][best_model]
                logger.info(f"The best predicted model by AlphaFold is {best_model} with a ptm of {iptm_ptm}")
            except Exception as e:
                logger.error("The scores are not available for the predicted models", exc_info=True)
                logger.info(e)
                
            file.close()
            """
            
            # Aphafold ranks automatically the models from best to worst
            best_model_file = f"ranked_0.pdb"
            logger.info(f"The file {best_model_file} will be used as a search model for molecular replacement!")

        # if there is an error with the json file the script will stop and give an error message, which will be logged
        except Exception as e:
            logger.error(f"The produced {infile} file in the AlphaFold file doesn't exist", exc_info=True)
            sys.exit(1)

        # AFpdbModel_path = infile = os.path.join(infile, f"{best_model_file}")
        # return the best model file
        return os.path.abspath(best_model_file)

    # adapted and changed version of phenix.process_predicted_model as below, which doesn't split the model into domains or chains
    # the formulas derived from the paper of Oeffner et al: https://journals.iucr.org/d/issues/2022/11/00/ai5009/ 
    def process_and_trim_pdb(self, output_dir):
        try:
            # read the pdb file and get the lines
            with open(f'{output_dir}/ranked_0.pdb', 'r') as file:
                lines = file.readlines()

            # loop over the lines and get the plddt score
            new_lines = []
            for line in lines:
                # Only process ATOM and HETATM lines
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    plddt = float(line[60:66].strip())
                    # Convert to scale of 0 to 1 if necessary
                    if plddt > 1:  
                        plddt /= 100
                    # Trimming residues with pLDDT < 70
                    if plddt < 0.7:  
                        continue
                    # Calculate delta or rmsd from the plddt score
                    delta = 1.5 * math.exp(4 * (0.7 - plddt)) 
                    # Calculate B-factor from delta
                    b_factor = (8 * math.pi**2 * delta**2) / 3
                    # Replace pseudo B-factor in line where plddt was
                    # The pseudo B-factor is a float with 6.2f format
                    # a new line is created with the new pseudo B-factor
                    new_line = line[:60] + f"{b_factor:6.2f}" + line[66:]
                    new_lines.append(new_line)

            # write the new lines to a new pdb file based on the new pseudo B-factor
            with open(f'{output_dir}/ranked_0_processed.pdb', 'w') as file:
                file.writelines(new_lines)
                
            shutil.move('job.log', output_dir)
            
            # if the processing of the model was successful the script will continue and log a message
            # inform the user that the processing of the model was successful and the next step will be molecular replacement (dimple)
            logger.info("Processing was succesful, starting molecular replacement...\n\n")
                
        # raise an error if the processing of the model was not successful
        except Exception as e:
            logger.error("Processing was not succesfull")
            logger.info(e)

    # method to process the best model from the prediction
    # the method process_predicted_model from phenix is used to process the model
    # the method takes the best model file and the output directory as arguments
    def phenix_process_predict(self, jobName, AFpdbModel_path, output_dir):
        script = "#!/usr/bin/env bash\n"

        # setting up the cluster environment or job specifications
        script += "#SBATCH --job-name=AF_processing\n"
        script += "#SBATCH --mem=0\n"
        script += "#SBATCH --time=01-00:00\n"

        # create a output and error file for the job
        script += "#SBATCH --output=phenix_%j.out\n"
        script += "#SBATCH --error=phenix_%j.err\n\n"

        # load the phenix module and purge the current modules
        script += "module purge\n"
        script += "module add gopresto Phenix/1.20.1-4487-Rosetta-3.10-PReSTO-8.0\n\n"

        # cd to the output directory where the processing of the model is done
        script += f"cd {output_dir}\n"

        # run phenix processing predicted model tool
        # Replace values in B-factor field with estimated B values.
        # Optionally remove low-confidence residues and split into domains.
        script += f"phenix.process_predicted_model {AFpdbModel_path} "

        # name the sbatch file and write the script to the file
        shellFile = f"{jobName}_slurm.sh"

        with open(shellFile,"w") as file:
            file.write(script)
            file.close()

        # get the job id from the monitor and print the status of the job to stdout
        job_id = UtilsMonitor.monitor_job(script=f"{jobName}_slurm.sh", name=f"Processing ({jobName})")

        # check if the processing of the model was successful
        # only if the job id is an integer the processing of the model was successful otherwise there is no job id
        try:
            int(job_id)
        except Exception as e:
            logger.error("Processing was not succesfull")
            logger.info(e)
   
        # move the slurm output files to the output directory
        move_file = [f"phenix_{job_id}.err", f"phenix_{job_id}.out"]
        for file in move_file:
            shutil.move(file, f"alf_output/{self.alphafold_job_id}")

        # if the processing of the model was successful the script will continue and log a message
        # inform the user that the processing of the model was successful and the next step will be molecular replacement (dimple)
        logger.info("Processing was succesful, starting molecular replacement...\n\n")

# for molecular replacement and refinement the program DIMPLE is used
# DIMPLE is a program for automated macromolecular crystal structure solution and takes in the processed model and the reflection data
# DIMPLE is available through: http://ccp4.github.io/dimple/ 
class mrALFO():
    """
    Molecular replacement with best predicted model from AlphaFold
    """
    def __init__(self, alphafold_job_id):
        self.alphafold_job_id = alphafold_job_id

    # DIMPLE - automated molecular replacement and refinement
    # the method takes the reflection data and the processed model as arguments
    def runDIMPLE(self, jobName, reflectionData_file, Model_file, output_dir):
        script = "#!/usr/bin/env bash\n"

        # setting up the cluster environment or job specifications
        script += "#SBATCH --job-name=AF_dimple\n"
        script += "#SBATCH --mem=0\n"
        script += "#SBATCH --time=01-00:00\n"

        # create a output and error file for the job
        script += "#SBATCH --output=dimple_%j.out\n"
        script += "#SBATCH --error=dimple_%j.err\n\n"

        # load the cpp4 module with dimple and purge the current modules
        script += "module purge\n"
        script += "module add gopresto CCP4/7.0.078-SHELX-ARP-8.0-17-PReSTO\n\n"

        # copy the reflection data to the output directory of dimple as it will be a subdirectory of the output directory of alphafold
        # cd to the output directory where the prediction is done
        script += f"cp {reflectionData_file} {output_dir}\n"
        script += f"cd {output_dir}\n\n"

        # run dimple with the processed model and reflection data and name the output directory dimpleMR
        script += f"dimple {Model_file} {reflectionData_file} dimpleMR"

        # name the sbatch file and write the script to the file
        shellFile = f"{jobName}_slurm.sh"

        with open(shellFile,"w") as file:
            file.write(script)
            file.close()

        # get the job id from the monitor and print the status of the job to stdout
        job_id = UtilsMonitor.monitor_job(script=f"{jobName}_slurm.sh", name=f"DIMPLE MR ({jobName})")

        # check if the prediction was successful
        # only if the job id is an integer the prediction was successful otherwise there is no job id
        try:
            int(job_id)
        except:
            logger.error("Molecular replacement with DIMPLE was not succesfull")

        # move the slurm output files to the output directory
        move_file = [f"dimple_{job_id}.err", f"dimple_{job_id}.out"]
        for file in move_file:
            shutil.move(file, f"alf_output{self.alphafold_job_id}")
            
        # potential file was created, which is identical to the screen.log file
        if os.path.exists("alf_outputNone"):
            os.remove("alf_outputNone")

        # if the prediction was successful the script will continue and log a message
        # inform that the pipeline is finished
        logger.info("Molecular replacement was successfull, all done!")

# main function to run the pipeline
if __name__ == "__main__":
    # check if the fasta file exists
    # if not the script will stop and give an error message, which will be logged
    fasta_file = args.fasta_path
    fasta_path = Path(fasta_file)
    if not fasta_path.is_file():
        logger.error("Path to fasta file does not exist")
        sys.exit("Path to fasta file does not exist")

    mtz_file = args.reflectionData_path
    mtzPath = Path(mtz_file)
    if not mtzPath.is_file():
        logger.error("Path to mtz file does not exist")
        sys.exit("Path to mtz file does not exist")

    # predict the model as a slurm job
    # instantiate the class and run the prediction
    # get the output directory and the fasta name for the next step by running the method run of the prediction class
    if not args.pdb_path:
        ALFOpred = ALFOpred()
        output_dir, fasta_name, job_id, preset = ALFOpred.run(fasta_path)
    else:
        output_dir = os.path.dirname(args.pdb_path)
        
          # check if it's a fasta file and get the name for the job
        try:
            with open(args.fasta_path, mode="r") as file:
                line = file.read()

                # reads the first line of the fasta file and check if it starts with a '>'
                # if not the script will stop and give an error message, which will be logged
                if not line.startswith('>'):
                    logger.error("The input is not a fasta file!")
                    sys.exit(1)

                # if the first line starts with a '>', it will be counted how many '>' or sequences are in the file
                # by that the script can determine if it's a monomer or multimer prediction and c hange the preset accordingly
                else:
                    if line.count('>') == 1:
                        preset = "monomer"
                    else:
                        preset = "multimer"
                    # get the name of the fasta file without the extension
                    line = line.strip()
                    fasta_name = line[1:5]
        
        # if there is an error with the fasta file the script will stop and give an error message, which will be logged
        # most likely the file doesn't exist or can't be opened
        except Exception as e:
            logger.error(f"{args} can not be open or does not exist!", exc_info=True)
            sys.exit(1)


    if preset == "monomer":
        # choose the best model and process it
        # instantiate the class and run the method choose_model and process_predict
        if not args.pdb_path:
            if not args.just_predict or args.predict_process:
                logger.info("Processing the best model from the prediction is running...")
                procALFO = procALFO(alphafold_job_id = ALFOpred.job_id)
                mrALFO = mrALFO(alphafold_job_id = ALFOpred.job_id)
                AFpdbModel_path = procALFO.get_model(output_dir)
                
                procALFO.process_and_trim_pdb(output_dir=output_dir)
                # procALFO.process_predict(jobName=f"{fasta_name}_proc", AFpdbModel_path=AFpdbModel_path, output_dir=output_dir)

            else:
                logger.info("--just-predict is given, only the prediction was done!")
        else:
            # Extracting the number from the file path
            match = re.search(r"/(\d+)_", args.pdb_path)
            job_id = match.group(1) if match else None
            procALFO = procALFO(alphafold_job_id = job_id)
            mrALFO = mrALFO(alphafold_job_id = job_id)
            AFpdbModel_path = args.pdb_path
            
            procALFO.process_and_trim_pdb(output_dir=output_dir)
            # procALFO.process_predict(jobName=f"{fasta_name}_proc", AFpdbModel_path=AFpdbModel_path, output_dir=output_dir)


        # use dimple to do molecular replacement
        # get the processed model and reflection data as arguments
        # as the processed model has a new name the path to the processed model is created by changing the path name of the best model
        # instantiate the class and run the method runDIMPLE
        if not args.just_predict or not args.predict_process:
            processedModel_file = f'{output_dir}/ranked_0_processed.pdb'
            mrALFO.runDIMPLE(jobName=f"{fasta_name}_dimp", reflectionData_file=mtz_file, Model_file=processedModel_file, output_dir=output_dir)

    else:
        if not args.pdb_path:
            if not args.just_predict or args.predict_process:
                logger.info("Processing the best model from the prediction is running...")
                procALFO = procALFO(alphafold_job_id = ALFOpred.job_id)
                mrALFO = mrALFO(alphafold_job_id = ALFOpred.job_id)
                AFpdbModel_path = procALFO.get_model(output_dir)
                procALFO.process_and_trim_pdb(output_dir=output_dir)
                # procALFO.process_predict(jobName=f"{fasta_name}_proc", AFpdbModel_path=AFpdbModel_path, output_dir=output_dir)
            else:
                logger.info("--just-predict is given, only the prediction was done!")
        else:
            # Extracting the number from the file path
            match = re.search(r"/(\d+)_", args.pdb_path)
            job_id = match.group(1) if match else None
            logger.info("Processing the best model from the prediction is running...")
            procALFO = procALFO(alphafold_job_id = job_id)
            mrALFO = mrALFO(alphafold_job_id = job_id)
            procALFO.process_and_trim_pdb(output_dir=output_dir)

        # use dimple to do molecular replacement
        # get the processed model and reflection data as arguments
        # as the processed model has a new name the path to the processed model is created by changing the path name of the best model
        # instantiate the class and run the method runDIMPLE
        if not args.just_predict or not args.predict_process:
            processedModel_file = f'{output_dir}/ranked_0_processed.pdb'
            mrALFO.runDIMPLE(jobName=f"{fasta_name}_dimp", reflectionData_file=mtz_file, Model_file=processedModel_file, output_dir=output_dir)
            
    # tidy up temporary files
    # if the option --tidy is given the temporary files will be deleted
    # this includes slurm scripts 
    if args.tidy:
        logger.info("Tidying up temporary files...")
        # remove the slurm scripts
        slurm_scripts = [f"AFpred_{fasta_name}_slurm.sh", f"{fasta_name}_proc_slurm.sh", f"{fasta_name}_dimp_slurm.sh",
                         '*.err', '*.out', '*.log']
        for script in slurm_scripts:
            if os.path.exists(script):
                os.remove(script)
            else:
                continue
            
        logger.info("Temporary files are removed!")