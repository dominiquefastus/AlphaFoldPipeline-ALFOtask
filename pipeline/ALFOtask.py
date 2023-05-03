#!/usr/bin/env python3

"""
This is a automated python based SLURM job submitting pipeline for alphafold prediction.

Author:     D. Fastus
"""

from pathlib import Path
import subprocess
import logging
import time
import sys
import os

from utils import UtilsMonitor

# build a logger for the SLURM script
logging.basicConfig(filename=f'job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)

class ALFOpred():
    """
    Initiates a AlphaFold prediction SLURM job
    """

    def __init__(self) -> None:
        pass

    # create sbatch file for server
    def sbatch_AFpred(self, jobName, fastaName, fastaPath, partition = "v100", mem = 0, time = "01-00:00"):
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
        script += "module add OpenSSL/1.0 fosscuda/2020b AlphaFold\n\n"

        # define alphafold database directory
        script += "export ALPHAFOLD_DATA_DIR=/sw/pkg/miv/mx/db/alphafold-2021b\n\n"

        # others
        script += "export CWD=`pwd`\n"
        script += f"mkdir --parents alf_output/$SLURM_JOBID\n\n"

        # slurmtmp
        script += f"cp {fastaPath} /local/slurmtmp.$SLURM_JOBID\n"
        script += f"cd /local/slurmtmp.$SLURM_JOBID\n\n"

        # run alphafold
        script += f"""alphafold \\
        --fasta_paths={fastaName+'.fasta'} \\
        --max_template_date=2020-05-14 \\
        --output_dir=$CWD/alf_output/$SLURM_JOBID \\
        --data_dir=$ALPHAFOLD_DATA_DIR"""
        
        shellFile = f"{jobName}_slurm.sh"

        with open(shellFile,"w") as file:
            file.write(script)
            file.close()

    # check if the output are complete and have the right format
    def check_out(self, out_dir):
        # List of files to check for
        files_to_check = ['ranked_0.pdb', 'relaxed_model_1.pdb', 'result_model_1.pkl', 'unrelaxed_model_1.pdb']

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
                line = file.readline()

                if not line.startswith('>'):
                    logger.error("The input is not a fasta file!")
                    sys.exit(1)
                else:
                    line = line.strip()
                    fasta_name = line[1:5]

        except Exception as e:
            logger.error(f"{args} can not be open or does not exist!", exc_info=True)
            sys.exit(1)
        
        # create the sbatch file
        ALFOpred.sbatch_AFpred(jobName=f"AFpred_{fasta_name}", fastaName=fasta_name, fastaPath=args)

        # inform the user about the status by getting the jobID
        CommandLine = f"sbatch AFpred_{fasta_name}_slurm.sh"
        job_id = subprocess.run(CommandLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        job_id = job_id.stdout.decode("ascii").rstrip().split()[-1]

        try:
            int(job_id)
        except:
            logger.error("Job was not succesfull")

        print(f"Job ID: {job_id}")
        time.sleep(2)
        UtilsMonitor.monitor_job(job_id)

        ALFOpred.check_out(f"alf_output/{fasta_name}")

class mrALFO():
    """
    Molecular replacement with best predicted model from AlphaFold
    """

    def __init__(self) -> None:
        pass

    def choose_model(self, indir):
        pass



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

    ALFOpred = ALFOpred()
    ALFOpred.run(fastapath)