#!/usr/bin/env python3

# import the necessary packages to log and interact with the file system
import logging
import os

logging.basicConfig(filename=f'job.log', filemode='w', format='%(asctime)s - %(message)s', level=logging.DEBUG)
logger = logging.getLogger(__name__)

# function to check if the output are complete and have the right format
# it simply checks if the files are present in the directory by looping over the files_to_check list
def loop_files(output_dir, files_to_check):

    # Check if all files are present in the directory
    for file in files_to_check:
        file_path = os.path.join(output_dir, file)

        if not os.path.isfile(file_path):
            logger.error(f"The files are not successfully generated, the {file} is missing...", exc_info=True)
            return False