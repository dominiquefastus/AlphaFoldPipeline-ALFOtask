#!/usr/bin/env python3

import subprocess
import logging
import time

START_JOB_CODES = ["PENDING","RUNNING","REQUEUED", "RESIZING","SUSPENDED"]
END_JOB_CODES = ["BOOT_FAIL","CANCELLED","COMPLETED","DEADLINE","FAILED","NODE_FAIL","OUT_OF_MEMORY","TIMEOUT","REVOKED","PREEMPTED"]

logging.basicConfig(filename=f'job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# make a continuos check if the job is finished on the server
def monitor_job(script, name):
    start_time = time.perf_counter()

    job = subprocess.check_output(['sbatch', script])
    time.sleep(2)

    jobnum = (job.rstrip().split()[-1]).decode('ascii')
    print(f"Starting Job: {jobnum} [{name}]")

    while True:
        print("{0}{1}{2}".format('waiting for ',name,' ...'),end='\r')
        job = subprocess.check_output(['sacct','-j',jobnum]).decode('ascii')
        result = map(lambda x: x in job,START_JOB_CODES)

        if not True in result:
            break
        time.sleep(1)


    end_time = time.perf_counter()
    total_time = end_time - start_time

    print(120*"=")
    # logger.info(f"Job {job_id} {squeue_status}")
    logger.info(f"Job {job} finished in {total_time} seconds")

    return jobnum
