#!/usr/bin/env python3

"""
This is a automated python based SLURM job submitting pipeline for openfold prediction.

Author:     D. Fastus
"""

from pathlib import Path
import argparse
import logging
import shutil
import json
import sys
import os

# import module to monitor the SLURM jobs
# import module to check the output

# sys.insert
from utils import UtilsMonitor

# build a logger for the SLURM script
logging.basicConfig(filename=f'alf_output/job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

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

class ALFOpred():
    """
    Initiates a AlphaFold prediction SLURM job
    """

    def __init__(self) -> None:
        pass

    # create sbatch file for server to run alphafold with monomer or multimer preset
    # uses the most stable and reliable version of alphafold on the server: AlphaFold/2.1.1
    def sbatch_AFpred(self, preset, jobName, fastaName, fastaPath, partition = "v100", mem = 0, time = "01-00:00"):
        script = "#!/usr/bin/env bash\n"

        # setting up the cluster environment or job specifications
        script += f"#SBATCH --job-name=AF_{fastaName}\n"

        if partition is not None:
            script += f"#SBATCH --partition={partition}\n"
        
        script += f"#SBATCH --mem={mem}\n"
        script += f"#SBATCH --time={time}\n"
        script += "#SBATCH --output=alphafold_%j.out\n"
        script += "#SBATCH --error=alphafold_%j.err\n\n"
        script += "#SBATCH --exclusive\n"


        # alphafold commands and requirements
        script += "module purge\n"
        script += "module add fosscuda/2020b AlphaFold\n\n"

        # define alphafold database directory
        script += "export ALPHAFOLD_DATA_DIR=/sw/pkg/miv/mx/db/alphafold-2021b\n\n"

        # others
        script += "export CWD=`pwd`\n"
        script += f"mkdir --parents alf_output/$SLURM_JOBID\n\n"

        # slurmtmp
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
        
        shellFile = f"{jobName}_slurm.sh"

        with open(shellFile,"w") as file:
            file.write(script)
            file.close()

    # check if the output are complete and have the right format
    def check_out(self, out_dir):
        # List of files to check for
        files_to_check = ['ranked_0.pdb', 'relaxed_model_1.pdb', 'result_model_1.pkl', 'unrelaxed_model_1.pdb', "ranking_debug.json"]

        # Check if all files are present in the directory
        for file in files_to_check:
            file_path = os.path.join(out_dir, file)

            if not os.path.isfile(file_path):
                logger.error(f"The files are not successfully generated, the {file} is missing...", exc_info=True)
                return False
            
        logger.info(f"The files are successfully generated in {out_dir}")
        return True

    def run(self, args):
        # check if it's a fasta file and get the name for the job
        try:
            with open(args,mode="r") as file:
                line = file.read()

                if not line.startswith('>'):
                    logger.error("The input is not a fasta file!")
                    sys.exit(1)
                else:
                    if line.count('>') == 1:
                        preset = "monomer"
                    else:
                        preset = "multimer"

                    line = line.strip()
                    fasta_name = line[1:5]

        except Exception as e:
            logger.error(f"{args} can not be open or does not exist!", exc_info=True)
            sys.exit(1)
        
        # create the sbatch file depenjding on the fasta file
        ALFOpred.sbatch_AFpred(preset=preset, jobName=f"AFpred_{fasta_name}", fastaName=fasta_name, fastaPath=args)
    

        # get the job id from the monitor and print the status of the job to stdout
        job_id = UtilsMonitor.monitor_job(script=f"AFpred_{fasta_name}_slurm.sh", name=f"Predicting ({fasta_name})")

        # move the slurm output files to the output directory
        # check if alphafold output is complete (files are present)
        move_file = [f"alphafold_{job_id}.err", f"alphafold_{job_id}.out"]
        for file in move_file:
            shutil.move(file, output_dir)


        output_dir = f"alf_output/{job_id}/{fasta_name}"
        ALFOpred.check_out(output_dir)

        logger.info(f"The prediction was successful, starting processing of model...\n\n")

        return output_dir, fasta_name


class procALFO():
    """
    AlphaFold predicted model selection and preprocessing
    """

    def __init__(self) -> None:
        pass

    def choose_model(self, infile):
        try:
            json_file = os.path.join(infile, "ranking_debug.json")
            file = open(json_file)
            ranking = json.load(file)

            best_model = ranking["order"][0]
            plddts = ranking["plddts"][best_model]
            best_model_file = f"ranked_0.pdb"

            file.close()

            logger.info(f"The best predicted model by AlphaFold is {best_model} with a plddts of {plddts}")
            logger.info(f"The file {best_model_file} will be used as a search model for molecular replacement!")

        except Exception as e:
            logger.error(f"The produced {infile} file in the AlphaFold file doesn't exist", exc_info=True)
            sys.exit(1)

        # AFpdbModel_path = infile = os.path.join(infile, f"{best_model_file}")
        return best_model_file

    def process_predict(self, jobName, AFpdbModel_path, output_dir):
        script = "#!/usr/bin/env bash\n"

        # setting up the cluster environment or job specifications
        script += "#SBATCH --job-name=AF_processing\n"
        script += "#SBATCH --mem=0\n"
        script += "#SBATCH --time=01-00:00\n"

        script += "#SBATCH --output=alphafold_%j.out\n"
        script += "#SBATCH --error=alphafold_%j.err\n\n"

         # alphafold commands and requirements
        script += "module purge\n"
        script += "module add gopresto Phenix/1.20.1-4487-Rosetta-3.10-PReSTO-8.0\n\n"

         # slurmtmp
        # script += f"cp {AFpdbModel_path} {output_dir}\n"
        script += f"cd {output_dir}\n"

        # run phenix processing predicted model tool
        # Replace values in B-factor field with estimated B values.
        # Optionally remove low-confidence residues and split into domains.
        script += f"phenix.process_predicted_model {AFpdbModel_path} "

        shellFile = f"{jobName}_slurm.sh"

        with open(shellFile,"w") as file:
            file.write(script)
            file.close()

        job_id = UtilsMonitor.monitor_job(script=f"{jobName}_slurm.sh", name=f"Processing ({jobName})")

        try:
            int(job_id)
        except:
            logger.error("Processing was not succesfull")

        move_file = [f"alphafold_{job_id}.err", f"alphafold_{job_id}.out"]
        for file in move_file:
            shutil.move(file, f"alf_output/{job_id}")

        logger.info("Processing was succesful, starting molecular replacement...\n\n")


class mrALFO():
    """
    Molecular replacement with best predicted model from AlphaFold
    """

    def __init__(self) -> None:
        pass

    # DIMPLE - automated refinement and ligand screening
    def runDIMPLE(self, jobName, reflectionData_file, processedModel_file, output_dir):
        script = "#!/usr/bin/env bash\n"

        # setting up the cluster environment or job specifications
        script += "#SBATCH --job-name=AF_processing\n"
        script += "#SBATCH --mem=0\n"
        script += "#SBATCH --time=01-00:00\n"

        script += "#SBATCH --output=alphafold_%j.out\n"
        script += "#SBATCH --error=alphafold_%j.err\n\n"

        # alphafold commands and requirements
        script += "module purge\n"
        script += "module add gopresto CCP4/7.0.078-SHELX-ARP-8.0-17-PReSTO\n\n"

            # slurmtmp
        script += f"cp {reflectionData_file} {output_dir}\n"
        script += f"cd {output_dir}\n\n"

        # run phenix processing predicted model tool
        # Replace values in B-factor field with estimated B values.
        # Optionally remove low-confidence residues and split into domains.
        script += f"dimple {processedModel_file} {reflectionData_file} dimpleMR"

        shellFile = f"{jobName}_slurm.sh"

        with open(shellFile,"w") as file:
            file.write(script)
            file.close()

        job_id = UtilsMonitor.monitor_job(script=f"{jobName}_slurm.sh", name=f"DIMPLE MR ({jobName})")

        try:
            int(job_id)
        except:
            logger.error("Molecular replacement with DIMPLE was not succesfull")

        move_file = ["job.log", f"alphafold_{job_id}.err", f"alphafold_{job_id}.out"]
        for file in move_file:
            shutil.move(file, f"alf_output/{job_id}")

        logger.info("Molecular replacement was successfull, all done!")

if __name__ == "__main__":
    # check for positional argument as fasta file
    if len(sys.argv) < 2:
        logger.error("Missing fasta file..")
        sys.exit("Error missing argument! Please provide a fasta file")
    args = sys.argv[1]
    fastapath = Path(args)
    if not fastapath.is_file():
        logger.error("Path to fasta file does not exist")
        sys.exit("Path to fasta file does not exist")

    elif len(sys.argv) < 3:
        logger.error("Missing mtz file..")
        sys.exit("Error missing argument! Please provide a mtz file")
    mtz_file = sys.argv[2]
    mtzPath = Path(mtz_file)
    if not mtzPath.is_file():
        logger.error("Path to mtz file does not exist")
        sys.exit("Path to mtz file does not exist")

    # predict the model as a slurm job
    ALFOpred = ALFOpred()
    output_dir, fasta_name = ALFOpred.run(fastapath)

    # choose the best model and process it
    procALFO = procALFO()
    AFpdbModel_path = procALFO.choose_model(output_dir)
    procALFO.process_predict(jobName=f"{fasta_name}_proc", AFpdbModel_path=AFpdbModel_path, output_dir=output_dir)

    # use dimple to do molecular replacement
    processedModel_file = AFpdbModel_path.replace(".pdb","_processed.pdb")
    mrALFO = mrALFO()
    mrALFO.runDIMPLE(jobName=f"{fasta_name}_dimp", reflectionData_file=mtz_file, processedModel_file=processedModel_file, output_dir=output_dir)
