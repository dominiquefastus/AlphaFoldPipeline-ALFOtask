"""
This is a automated python based SLURM job submitting pipeline for alphafold prediction.

Author:     D. Fastus
"""

import subprocess
import time
import sys
import os

# create sbatch file for server
def sbatch_AFpred(jobName, fastaName, fastaPath, partition = "v100", mem = 0, time = "01-00:00"):
    script = "#!/bin/bash\n"

    # setting up the cluster environment or job specifications
    script += f"#SBATCH --job-name=AF_{fastaName}\n"

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
    script += f"cp {fastaPath} /local/slurmtmp.$SLURM_JOBID\n"
    script += f"cd /local/slurmtmp.$SLURM_JOBID\n\n"

    # run alphafold
    script += f"""alphafold \\
                --fasta_paths={fastaName+'.fasta'} \\
                --max_template_date=2020-05-14 \\
                --output_dir=$CWD/alf_output \\
                --data_dir=$ALPHAFOLD_DATA_DIR
               """
    
    shellFile = f"{jobName}_slurm.sh"

    with open(shellFile,"w") as file:
        file.write(script)
        file.close()

# make a continuos check if the job is finished on the server
def monitor_job(job_id):
    start_time = time.time()
    while True:
        # Run the squeue command and capture the output
        squeue_output = subprocess.run(f"sacct -j {job_id}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        squeue_output = squeue_output.stdout.decode("ascii").rstrip().split()

        print("Job is still running...")
        time.sleep(20)  # Wait for 5 minutes before checking again

        # Check if the job ID is still present in the output
        if "FAILED" in squeue_output:
            print(f"The job {job_id} failed to run")
            break
        elif "COMPLETED" in squeue_output:
            print(f"The job {job_id} is completed")
            break

    end_time = time.time()
    total_time = end_time - start_time

    with open('job.log', 'a') as f:
        f.write(f"Job {job_id} finished in {total_time} seconds\n")

    print(f"Job {job_id} has finished in {total_time} seconds")


# check if the output are complete and have the right format
def check_out(out_dir):
    # List of files to check for
    files_to_check = ['ranked_0.pdb', 'relaxed_model_1.pdb', 'result_model_1.pkl', 'unrelaxed_model_1.pdb']

    # Check if all files are present in the directory
    for file in files_to_check:
        file_path = os.path.join(out_dir, file)

        if not os.path.isfile(file_path):
            print(f"The files are not successfully generated, the {file} is missing...")
            return False
        
    print(f"The files are successfully generated in {out_dir}")
    return True

def main(args):
    if len(args) < 2 or len(args) == 0:
        print("Please provide only one fasta file")
        sys.exit(1)

    # che if it's a fasta file and get the name for the job
    try:
        with open(args,mode="r") as file:
            line = file.readline()

            if not line.startswith('>'):
                print("The input is not a fasta file!")
                sys.exit(1)
            else:
                line = line.strip()
                fasta_name = line[1:].split("|")[0]

    except FileNotFoundError:
        print(f"{args} can not be open or does not exist!")
        sys.exit(1)
    
    # create the sbatch file
    sbatch_AFpred(jobName=f"AFpred_{fasta_name}", fastaName=fasta_name, fastaPath=os.getcwd()+args)

    # inform the user about the status by getting the jobID
    CommandLine = f"sbatch --parsable AFpred_{fasta_name}_slurm.sh"
    job_id = subprocess.run(CommandLine, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    job_id = job_id.stdout.decode("ascii").rstrip().split()[-1]

    try:
        int(job_id)
    except:
        print("Job was not succesfull")

    monitor_job(job_id)

    check_out(f"alf_output/{fasta_name}")

if __name__ == "__main__":
    main(sys.argv[1])