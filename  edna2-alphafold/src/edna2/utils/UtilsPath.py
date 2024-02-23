#
# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__authors__ = ["O. Svensson"]
__license__ = "MIT"
__date__ = "21/04/2019"

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# mxv1/src/EDHandlerESRFPyarchv1_0.py

import os
import time
import pathlib
import tempfile
import shutil
import hashlib
import string
import random

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

logger = UtilsLogging.getLogger()

DEFAULT_TIMEOUT = 120  # s
STRING_CHARACTERS = string.ascii_lowercase + string.digits + "_"
STRING_LENGTH = 8


def getWorkingDirectory(task, inData, workingDirectorySuffix=None):
    parentDirectory = inData.get("workingDirectory", None)
    if parentDirectory is None:
        parentDirectory = os.getcwd()
    parentDirectory = pathlib.Path(parentDirectory)
    if workingDirectorySuffix is None:
        # Create unique directory
        # workingDirectory = tempfile.mkdtemp(
        #     prefix=task.__class__.__name__ + "_", dir=parentDirectory
        # )

        # os.chmod(workingDirectory, 0o755)
        # workingDirectory = pathlib.Path(workingDirectory)
        workingDirectory = makeRandomDirectoryPath(
            prefix=task.__class__.__name__, dir=parentDirectory
        )
        workingDirectory.mkdir(mode=0o775, parents=True, exist_ok=False)

    else:
        # Here we assume that the user knows what he is doing and there's no
        # race condition for creating the working directory!
        workingDirectoryName = (
            task.__class__.__name__ + "_" + str(workingDirectorySuffix)
        )
        workingDirectory = parentDirectory / workingDirectoryName
        index = 1
        while workingDirectory.exists():
            workingDirectoryName = (
                task.__class__.__name__
                + "_"
                + str(workingDirectorySuffix)
                + "_{0:02d}".format(index)
            )
            workingDirectory = parentDirectory / workingDirectoryName
            index += 1
        workingDirectory.mkdir(mode=0o775, parents=True, exist_ok=False)
    workingDirectory = stripDataDirectoryPrefix(workingDirectory)
    return workingDirectory


def makeRandomDirectoryPath(prefix, dir):
    "returns a random directory path that does not exist yet"
    randomString = "".join(
        random.choice(STRING_CHARACTERS) for _ in range(STRING_LENGTH)
    )
    dirPath = prefix + "_" + randomString
    returnPath = pathlib.Path(dir) / dirPath
    while returnPath.exists():
        randomString = "".join(
            random.choice(STRING_CHARACTERS) for _ in range(STRING_LENGTH)
        )
        dirPath = prefix + "_" + randomString
        returnPath = pathlib.Path(dir) / dirPath
    return returnPath


def createPyarchFilePath(filePath):
    """
    This method translates from an ESRF "visitor" path to a "pyarch" path:
    /data/visitor/mx415/id14eh1/20100209 -> /data/pyarch/2010/id14eh1/mx415/20100209
    """
    pyarchFilePath = None
    if isinstance(filePath, str):
        filePath = pathlib.Path(filePath)
    listOfDirectories = list(filePath.parts)
    if UtilsConfig.isMAXIV():
        # NOTE: Temporary fix that should already be in the branch. This is just a temporary
        # patch to avoid checking out the branch in the current state.
        if "visitors" in filePath.parts:
            index = filePath.parts.index("visitors")
        elif "proprietary" in filePath.parts:
            index = filePath.parts.index("proprietary")
        else:
            logger.error(
                "The /data/ directory should contain either visitors/ or proprietary/\n"
            )
            assert (
                False
            ), "Unexpected filesystem dirs"  # Don't continue if filesystem is not setup right.
        return pathlib.Path("/data/staff/ispybstorage").joinpath(
            *filePath.parts[index:]
        )

    if UtilsConfig.isEMBL():
        if "p13" in listOfDirectories[0:3] or "P13" in listOfDirectories[0:3]:
            pyarchFilePath = os.path.join("/data/ispyb/p13", *listOfDirectories[4:])
        else:
            pyarchFilePath = os.path.join("/data/ispyb/p14", *listOfDirectories[4:])
        return pyarchFilePath
    listBeamlines = [
        "bm07",
        "id23eh1",
        "id23eh2",
        "id29",
        "id30a1",
        "id30a2",
        "id30a3",
        "id30b",
    ]

    if (
        "data" in listOfDirectories
        and len(listOfDirectories) > 5
        and listOfDirectories[1] != "data"
    ):
        while listOfDirectories[1] != "data" and len(listOfDirectories) > 5:
            del listOfDirectories[1]

    # Check that we have at least four levels of directories:
    if len(listOfDirectories) > 5:
        dataDirectory = listOfDirectories[1]
        secondDirectory = listOfDirectories[2]
        thirdDirectory = listOfDirectories[3]
        fourthDirectory = listOfDirectories[4]
        fifthDirectory = listOfDirectories[5]
        year = fifthDirectory[0:4]
        proposal = None
        beamline = None
        listOfRemainingDirectories = listOfDirectories
        if dataDirectory == "data" and secondDirectory == "gz":
            if thirdDirectory == "visitor":
                proposal = fourthDirectory
                beamline = fifthDirectory
            elif fourthDirectory == "inhouse":
                proposal = fifthDirectory
                beamline = thirdDirectory
            else:
                raise RuntimeError(
                    "Illegal path for UtilsPath.createPyarchFilePath: "
                    + "{0}".format(filePath)
                )
            listOfRemainingDirectories = listOfDirectories[6:]
        elif dataDirectory == "data" and secondDirectory == "visitor":
            proposal = listOfDirectories[3]
            beamline = listOfDirectories[4]
            listOfRemainingDirectories = listOfDirectories[5:]
        elif dataDirectory == "data" and secondDirectory in listBeamlines:
            beamline = secondDirectory
            proposal = listOfDirectories[4]
            listOfRemainingDirectories = listOfDirectories[5:]
        if proposal is not None and beamline is not None:
            pyarchFilePath = pathlib.Path("/data/pyarch") / year / beamline
            pyarchFilePath = pyarchFilePath / proposal
            for directory in listOfRemainingDirectories:
                pyarchFilePath = pyarchFilePath / directory
    if pyarchFilePath is None:
        logger.warning(
            "UtilsPath.createPyarchFilePath: path not converted for"
            + " pyarch: %s " % filePath
        )
    else:
        pyarchFilePath = pyarchFilePath.as_posix()
    return pyarchFilePath


def waitForFile(file, expectedSize=None, timeOut=DEFAULT_TIMEOUT):
    """Wait for the file to appear on disk."""
    file_path = pathlib.Path(file)
    file_size = None
    final_size = None
    has_timed_out = False
    should_continue = True
    fileHash = "abc"
    fileHash_old = "def"
    file_dir = file_path.parent
    if os.name != "nt" and file_dir.exists():
        # Patch provided by Sebastien 2018/02/09 for forcing NFS cache:
        # logger.debug("NFS cache clear, doing os.fstat on directory {0}".format(fileDir))
        fd = os.open(file_dir.as_posix(), os.O_DIRECTORY)
        stat_result = os.fstat(fd)
        os.close(fd)
        # logger.debug("Results of os.fstat: {0}".format(statResult))
    # Check if file is there
    if file_path.exists():
        file_size = file_path.stat().st_size
        file_mtime = file_path.stat().st_mtime
        time.sleep(0.1)
        # if expectedSize is not None:
        #     # Check size
        #     if file_size > expectedSize:
        #         should_continue = False
        # final_size = file_size
    else:
        file_size = 0
        file_mtime = 0
        time.sleep(0.1)
    if should_continue:
        logger.info("Waiting for file %s" % file_path)
        #
        time_start = time.time()
        while should_continue and not has_timed_out:
            if os.name != "nt" and file_dir.exists():
                # Patch provided by Sebastien 2018/02/09 for forcing NFS cache:
                # logger.debug("NFS cache clear, doing os.fstat on directory {0}".format(fileDir))
                fd = os.open(file_dir.as_posix(), os.O_DIRECTORY)
                stat_result = os.fstat(fd)  # noqa F841
                os.close(fd)
                # logger.debug("Results of os.fstat: {0}".format(statResult))
            time_elapsed = time.time() - time_start
            # Check if time out
            if time_elapsed > timeOut:
                has_timed_out = True
                str_warning = f"Timeout while waiting for file {file_path}"
                logger.warning(str_warning)
            else:
                # Check if file is there
                if file_path.exists():
                    file_size_new = file_path.stat().st_size
                    file_mtime_new = file_path.stat().st_mtime
                    if expectedSize is not None:
                        # Check that it has right size
                        if (
                            file_size > expectedSize
                            and file_size_new == file_size
                            and file_mtime_new == file_mtime
                        ):
                            should_continue = False
                    else:
                        if file_size_new == file_size and file_mtime_new == file_mtime:
                            should_continue = False
                    final_size = file_size
                    file_size = file_size_new
                    file_mtime == file_mtime_new
            if should_continue:
                # Sleep 1 s
                time.sleep(1)
    return has_timed_out, final_size


def get_md5Hash(file):
    file_path = pathlib.Path(file)
    if not file_path.exists():
        return None
    with open(file_path, "rb") as fp:
        md5Hash = hashlib.md5()
        chunk = fp.read(65536)
        while len(chunk) > 0:
            md5Hash.update(chunk)
            chunk = fp.read(65536)
    hexHash = md5Hash.hexdigest()
    return hexHash


def stripDataDirectoryPrefix(data_directory):
    """Removes any paths before /data/..., e.g. /gpfs/easy/data/..."""
    list_paths = str(data_directory).split(os.sep)
    if "data" in list_paths:
        while list_paths[1] != "data":
            list_paths = [list_paths[0]] + list_paths[2:]
        new_data_directory = os.sep.join(list_paths)
    else:
        new_data_directory = data_directory
    return pathlib.Path(new_data_directory)


def systemCopyFile(fp_in, fp_out):
    """Uses shutil.copy2 to copy files."""
    try:
        logger.debug(f"Copying {fp_in} to {fp_out}...")
        fout = shutil.copy2(fp_in, fp_out)
    except Exception as e:
        logger.error(f"Copying {fp_in} to {fp_out} failed: {e}.")
        fout = None
    return fout
