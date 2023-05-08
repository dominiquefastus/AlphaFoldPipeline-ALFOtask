#!/usr/bin/env python3

import logging
import os

START_JOB_CODES = ["PENDING","RUNNING","REQUEUED", "RESIZING","SUSPENDED"]
END_JOB_CODES = ["BOOT_FAIL","CANCELLED","COMPLETED","DEADLINE","FAILED","NODE_FAIL","OUT_OF_MEMORY","TIMEOUT","REVOKED","PREEMPTED"]

logging.basicConfig(filename=f'job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# check if the output are complete and have the right format
def check_out(self, out_dir, files_to_check):

    # Check if all files are present in the directory
    for file in files_to_check:
        file_path = os.path.join(out_dir, file)

        if not os.path.isfile(file_path):
            logger.error(f"The files are not successfully generated, the {file} is missing...", exc_info=True)
            return False
        
    logger.info(f"The files are successfully generated in {out_dir}")
    return True