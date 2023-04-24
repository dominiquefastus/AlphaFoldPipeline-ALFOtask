"""
This is a automated python based SLURM job submitting pipeline for alphafold prediction.
"""

import subprocess
import time
import glob
import os


# scan current working directory for fasta file to be used for AFpred
fasta_files = glob.glob(os.path.join(os.getcwd(), '*.fasta'))
if fasta_files:
    fasta_path = fasta_files[0]
    fasta_name = os.path.splitext(os.path.basename(fasta_path))[0]
    print(f"Found a fasta file: {fasta_path}")
else:
    print("No fasta files found in the current working directory.")
    exit()
    

# optional parse the fasta file for protein name, check structure
# possibility to estimate the prediction time from the residues
def protParser(fasta_file): 
    try:
        with open(fasta_file,mode="r") as file:
            for line in file:
                if line.startswith('>'):
                    line = line.strip()
                    protein_name = line[1:].split()[0]
                else:
                    num_residue = 0
                    for residue in line:
                      num_residue += 1
    except:
        pass
            
    return fasta_name, num_residue

# create sbatch file for server
def sbatch_AFpred(jobName = f"AFpred_{fasta_name}", partition = "v100", mem = 0, time = "01-00:00"):
    script = "#!/bin/bash\n"

    # setting up the cluster environment or job specifications
    script += f"#SBATCH --job-name=AF_{fasta_name}\n"

    if partition is not None:
        script += f"#SBATCH --partition={partition}\n"
    
    script += f"#SBATCH --exclusive\n"
    script += f"#SBATCH --mem={mem}\n"
    script += f"#SBATCH --time={time}\n"
    script += f"#SBATCH --output=alphafold_%j.out\n"
    script += f"#SBATCH --error=alphafold_%j.err\n\n"

    # alphafold commands and requirements
    script += "module purge\n"
    script += "module add fosscuda/2020b AlphaFold\n\n"

    # define alphafold database directory
    script += "export ALPHAFOLD_DATA_DIR=/sw/pkg/miv/mx/db/alphafold-2021b\n\n"

    # others
    script += "export CWD=`pwd`\n"
    script += f"mkdir --parents {jobName+'_output'}\n\n"

    # slurmtmp
    script += f"cp {fasta_path} /local/slurmtmp.$SLURM_JOBID\n"
    script += f"cd /local/slurmtmp.$SLURM_JOBID\n\n"

    # run alphafold
    script += f"""alphafold \\
                --fasta_paths={fasta_name+'.fasta'} \\
                --max_template_date=2020-05-14 \\
                --output_dir=$CWD/alf_output \\
                --data_dir=$ALPHAFOLD_DATA_DIR
               """
    
    shellFile = f"{jobName}_slurm.sh"

    with open(shellFile,"w") as file:
        file.write(script)
        file.close()
    
    # make the sbatch file executable
    st = os.stat(shellFile)
    os.chmod(shellFile, st.st_mode | 0o755)

# make a continuos check if the job is finished on the server
def monitor_job(job_id):
    start_time = time.time()
    while True:
        # Run the squeue command and capture the output
        squeue_output = subprocess.run(f"sacct -j {job_id}", shell=True)

        # Check if the job ID is still present in the output
        if "RUNNING" or "WAITING" in squeue_output:
            print("Job is still running...")
            time.sleep(300)  # Wait for 5 minutes before checking again

        else:
            end_time = time.time()
            total_time = end_time - start_time

            with open('job.log', 'a') as f:
                f.write(f"Job {job_id} finished in {total_time} seconds\n")

            print(f"Job {job_id} has finished in {total_time} seconds")
            break

# since the sbatch file is already running this function is unneccessary, need to be checked
def run(shellFile):
    SlurmCommandLine = f"sbatch {shellFile}"

    # execute the sbatch script
    subprocess.run(SlurmCommandLine, shell=True)

# check if the output are complete and have the right format
def check_out(out_dir):
    
    # List of files to check for
    files_to_check = ['ranked_0.pdb', 'relaxed_model_1.pdb', 'result_model_1.pkl', 'unrelaxed_model_1.pdb']

    # Check if all files are present in the directory
    for file in files_to_check:
        file_path = os.path.join(out_dir, file)

        if not os.path.isfile(file_path):
            os.makedirs("empty_folder")
            print(f"The files are not successfully generated, the {file} is missing...")
            return False
        
    print(f"The files are successfully generated in {out_dir}")
    return True

# get the name for the job name 
fasta_name, num_residues = protParser(fasta_path)

# create the sbatch file
sbatch_AFpred()

# run the sbatch file or SLURM JOB
# run(f"AFpred_{fasta_name}_slurm.sh")

# inform the user about the status
CommandLine = f"sbatch --parsable AFpred_{fasta_name}_slurm.sh"
job_id = subprocess.run(CommandLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
monitor_job(job_id)

check_out(f"alf_output/7QRZ")