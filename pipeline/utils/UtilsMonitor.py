#!/usr/bin/env python3

import subprocess
import logging
import time

START_JOB_CODES = ["PENDING","RUNNING","REQUEUED", "RESIZING","SUSPENDED"]
END_JOB_CODES = ["BOOT_FAIL","CANCELLED","COMPLETED","DEADLINE","FAILED","NODE_FAIL","OUT_OF_MEMORY","TIMEOUT","REVOKED","PREEMPTED"]

logging.basicConfig(filename=f'job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# make a continuos check if the job is finished on the server
def monitor_job(job_id, name):
    start_time = time.perf_counter()


# Set the default status to running
    squeue_status = "RUNNING"
    time.sleep(10)

    while squeue_status in START_JOB_CODES:
        # Run the squeue command and capture the output
        squeue_output = subprocess.run(f"sacct -j {job_id} -o State", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        squeue_output = squeue_output.stdout.decode("ascii").rstrip().split()
        squeue_status = squeue_output[2]
        print(squeue_status)

        print(f"Job {job_id} [{name}] is still running...")
        time.sleep(10)  # Wait for 5 minutes before checking again

        # Check if the job ID is still present in the output

    end_time = time.perf_counter()
    total_time = end_time - start_time

    print(120*"=")
    logger.info(f"Job {job_id} {squeue_status}")
    logger.info(f"Job {job_id} finished in {total_time} seconds")