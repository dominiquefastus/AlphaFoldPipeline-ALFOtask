#!/usr/bin/env python3

# import the necessary packages to log, measure time and run subprocesses
import subprocess
import logging
import time

# define slurm job codes to check if the job is running or finished
START_JOB_CODES = ["PENDING","RUNNING","REQUEUED", "RESIZING","SUSPENDED"]
# END_JOB_CODES = ["BOOT_FAIL","CANCELLED","COMPLETED","DEADLINE","FAILED","NODE_FAIL","OUT_OF_MEMORY","TIMEOUT","REVOKED","PREEMPTED"]

# set up the logger
logging.basicConfig(filename=f'job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# make a continuos check if the job is finished on the server
# the function takes as input the script to run and the name of the job
def monitor_job(script, name):
    # start the timer to measure the time needed to run the job
    start_time = time.perf_counter()

    # submit the job to the server with the sbatch command
    job = subprocess.check_output(['sbatch', script])
    # wait 2 seconds to let the job start
    time.sleep(2)

    # get the job id from the output of the sbatch command
    # the job id is the last element of the output
    # the output is a byte string, so we need to decode it to ascii
    job_id = (job.rstrip().split()[-1]).decode('ascii')
    # print the job id and the name of the job
    print(f"Starting Job: {job_id} [{name}]")

    # check the status of the job on the server and break the loop when the job is finished
    # wait 1 second between each check
    while True:
        # print a waiting phrase to let the user know that the job is running
        print("waiting for {0}{1}".format(name,' ...'), end='\r')
        # display accoutning information for the job and decode the output to ascii
        job = subprocess.check_output(['sacct','-j',job_id]).decode('ascii')
        # map the job status to over the job codes to check if the job is finished
        # applies the lambda function to each element of the codes list and returns a list of booleans
        # indicat the presence of the job code in the job status
        result = map(lambda x: x in job,START_JOB_CODES)

        # if the job is finished, break the loop
        if not True in result:
            break
        # wait 1 second between each check
        time.sleep(1)

    # end the timer to measure the time needed to run the job
    # calculate the total time needed to run the job
    end_time = time.perf_counter()
    total_time = end_time - start_time

    print(120*"=")
    # log the job id and the name of the job and the time needed to run the job
    logger.info(f"Job {job} finished in {total_time} seconds")

    # return the job id for further use
    return job_id
