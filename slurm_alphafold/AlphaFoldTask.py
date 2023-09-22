#!/usr/bin/env python3

"""
This is a automated python based SLURM job submitting pipeline for openfold prediction.

Author:     D. Fastus
"""

# necessary imports for the script
from pathlib import Path
import argparse
import logging
import shutil
import json
import sys
import os

# import module to monitor the SLURM jobs
from utils import UtilsMonitor

# build a logger for the SLURM script
logging.basicConfig(filename=f'alf_output/job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# add a handler to print the log messages to stdout
# also define the format of the log messages
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

'''
# Add argparse to the script to get diffent input files
parser = argparse.ArgumentParser(description="AlphaFold prediction pipeline")

# Add arguments
parser.add_argument("-f", "--fasta", help="Fasta file to predict protein structure")
parser.add_argument("-p", "--pdb", help="Pdb file to extract sequence and predict protein structure")
parser.add_argument("-m", "--mtz", help="Mtz file to do molecular replacement")
parser.add_argument("-o", "--output", help="Output directory for the results", default="alf_output", required=False)

# Parse arguments  
args = parser.parse_args()
'''

# the pipeline is designed in a object oriented way to keep the code clean and readable but als keep the structure 
# of the edna framework for easier implementation 

# the pipeline is divided in three steps represented by three classes
# class ALFOpred: initiate the alphafold prediction as a slurm job by using the installed alphafold module on the server
# alphafold is a deep learning based method for protein structure prediction and available through: 
# https://github.com/google-deepmind/alphafold 
class ALFOpred():
    """
    Initiates a AlphaFold prediction SLURM job
    """

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

        # create aa output and error file for the job
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
    def check_out(self, out_dir):
        # List of files to check for
        files_to_check = ['ranked_0.pdb', 'relaxed_model_1.pdb', 'result_model_1.pkl', 'unrelaxed_model_1.pdb', "ranking_debug.json"]

        # Check if all defined files are present in the directory
        for file in files_to_check:
            file_path = os.path.join(out_dir, file)

            # if the file is not present in the directory the script will stop and give an error message, which will be logged
            if not os.path.isfile(file_path):
                logger.error(f"The files are not successfully generated, the {file} is missing...", exc_info=True)
                return False

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

        # get the output directory to check if the prediction was successful and move some files
        output_dir = f"alf_output/{job_id}/{fasta_name}"

        # move the slurm output files to the output directory
        move_file = [f"alphafold_{job_id}.err", f"alphafold_{job_id}.out"]
        for file in move_file:
            shutil.move(file, output_dir)

        # check if the prediction was successful
        # check if alphafold output is complete (files are present)
        ALFOpred.check_out(output_dir)

        # if the prediction was successful the script will continue and log a message
        logger.info(f"The prediction was successful, starting processing of model...\n\n")

        # return the output directory and the fasta name for the next step
        return output_dir, fasta_name


# the model from the prediction is processed and the best model is selected
# for further molecular replacement and refinement the model is processed with phenix
# phenix is a software for the automated determination of molecular structures using X-ray crystallography and other methods
# the method process_predicted_model is used to process the model and is available through: https://www.phenix-online.org/documentation/reference/process_predicted_model.html
# basically the method replaces the B-factor values with estimated B values and removes low-confidence residues and split into domains if needed
class procALFO():
    """
    AlphaFold predicted model selection and preprocessing
    """

    # method to choose the best model from the prediction
    # the best model is chosen by the plddts score
    # this is only necessary for loging purposes and to give the user information about the best model
    # alphafold ranks the models automatically from best to worst and the best model is always ranked as 0
    def choose_model(self, infile):
        try:
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
            else:
                logger.error("The scores are not available for the predicted models", exc_info=True)

            # Aphafold ranks automatically the models from best to worst
            best_model_file = f"ranked_0.pdb"
            logger.info(f"The file {best_model_file} will be used as a search model for molecular replacement!")

            file.close()

        # if there is an error with the json file the script will stop and give an error message, which will be logged
        except Exception as e:
            logger.error(f"The produced {infile} file in the AlphaFold file doesn't exist", exc_info=True)
            sys.exit(1)

        # AFpdbModel_path = infile = os.path.join(infile, f"{best_model_file}")
        # return the best model file
        return best_model_file

    # method to process the best model from the prediction
    # the method process_predicted_model from phenix is used to process the model
    # the method takes the best model file and the output directory as arguments
    def process_predict(self, jobName, AFpdbModel_path, output_dir):
        script = "#!/usr/bin/env bash\n"

        # setting up the cluster environment or job specifications
        script += "#SBATCH --job-name=AF_processing\n"
        script += "#SBATCH --mem=0\n"
        script += "#SBATCH --time=01-00:00\n"

        # create a output and error file for the job
        script += "#SBATCH --output=alphafold_%j.out\n"
        script += "#SBATCH --error=alphafold_%j.err\n\n"

        # load the phenix module and purge the current modules
        script += "module purge\n"
        script += "module add gopresto Phenix/1.20.1-4487-Rosetta-3.10-PReSTO-8.0\n\n"

        # cd to the output directory where the prediction is done
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

        # check if the prediction was successful
        # only if the job id is an integer the prediction was successful otherwise there is no job id
        try:
            int(job_id)
        except:
            logger.error("Processing was not succesfull")

        # move the slurm output files to the output directory
        move_file = [f"alphafold_{job_id}.err", f"alphafold_{job_id}.out"]
        for file in move_file:
            shutil.move(file, f"alf_output/{job_id}")

        # if the prediction was successful the script will continue and log a message
        # inform the user that the prediction was successful and the next step will be molecular replacement (dimple)
        logger.info("Processing was succesful, starting molecular replacement...\n\n")

# for molecular replacement and refinement the program DIMPLE is used
# DIMPLE is a program for automated macromolecular crystal structure solution and takes in the processed model and the reflection data
# DIMPLE is available through: http://ccp4.github.io/dimple/ 
class mrALFO():
    """
    Molecular replacement with best predicted model from AlphaFold
    """

    # DIMPLE - automated molecular replacement and refinement
    # the method takes the reflection data and the processed model as arguments
    def runDIMPLE(self, jobName, reflectionData_file, processedModel_file, output_dir):
        script = "#!/usr/bin/env bash\n"

        # setting up the cluster environment or job specifications
        script += "#SBATCH --job-name=AF_dimple\n"
        script += "#SBATCH --mem=0\n"
        script += "#SBATCH --time=01-00:00\n"

        # create a output and error file for the job
        script += "#SBATCH --output=alphafold_%j.out\n"
        script += "#SBATCH --error=alphafold_%j.err\n\n"

        # load the cpp4 module with dimple and purge the current modules
        script += "module purge\n"
        script += "module add gopresto CCP4/7.0.078-SHELX-ARP-8.0-17-PReSTO\n\n"

        # copy the reflection data to the output directory of dimple as it will be a subdirectory of the output directory of alphafold
        # cd to the output directory where the prediction is done
        script += f"cp {reflectionData_file} {output_dir}\n"
        script += f"cd {output_dir}\n\n"

        # run dimple with the processed model and reflection data and name the output directory dimpleMR
        script += f"dimple {processedModel_file} {reflectionData_file} dimpleMR"

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

        # move the slurm output files to the output directory and also the log file 
        move_file = ["job.log", f"alphafold_{job_id}.err", f"alphafold_{job_id}.out"]
        for file in move_file:
            shutil.move(file, f"alf_output/{job_id}")

        # if the prediction was successful the script will continue and log a message
        # inform that the pipeline is finished
        logger.info("Molecular replacement was successfull, all done!")

# main function to run the pipeline
if __name__ == "__main__":
    # check for positional argument as fasta file
    # check if the fasta file exists
    # if not the script will stop and give an error message, which will be logged
    if len(sys.argv) < 2:
        logger.error("Missing fasta file..")
        sys.exit("Error missing argument! Please provide a fasta file")

    args = sys.argv[1]
    fastapath = Path(args)
    if not fastapath.is_file():
        logger.error("Path to fasta file does not exist")
        sys.exit("Path to fasta file does not exist")
    
    # check for positional argument as mtz file
    # check if the mtz file exists
    # if not the script will stop and give an error message, which will be logged
    elif len(sys.argv) < 3:
        logger.error("Missing mtz file..")
        sys.exit("Error missing argument! Please provide a mtz file")
    mtz_file = sys.argv[2]
    mtzPath = Path(mtz_file)
    if not mtzPath.is_file():
        logger.error("Path to mtz file does not exist")
        sys.exit("Path to mtz file does not exist")

    # predict the model as a slurm job
    # instantiate the class and run the prediction
    # get the output directory and the fasta name for the next step by running the method run of the prediction class
    ALFOpred = ALFOpred()
    output_dir, fasta_name = ALFOpred.run(fastapath)

    # choose the best model and process it
    # instantiate the class and run the method choose_model and process_predict
    procALFO = procALFO()
    AFpdbModel_path = procALFO.choose_model(output_dir)
    procALFO.process_predict(jobName=f"{fasta_name}_proc", AFpdbModel_path=AFpdbModel_path, output_dir=output_dir)

    # use dimple to do molecular replacement
    # get the processed model and reflection data as arguments
    # as the processed model has a new name the path to the processed model is created by changing the path name of the best model
    # instantiate the class and run the method runDIMPLE
    processedModel_file = AFpdbModel_path.replace(".pdb","_processed.pdb")
    mrALFO = mrALFO()
    mrALFO.runDIMPLE(jobName=f"{fasta_name}_dimp", reflectionData_file=mtz_file, processedModel_file=processedModel_file, output_dir=output_dir)
