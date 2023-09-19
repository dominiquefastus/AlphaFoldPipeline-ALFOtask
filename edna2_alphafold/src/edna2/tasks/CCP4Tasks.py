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
__date__ = "10/05/2019"

# Corresponding EDNA code:
# https://github.com/olofsvensson/edna-mx
# mxPluginExec/plugins/EDPluginGroupCCP4-v1.0/plugins/EDPluginExecAimlessv1_0.py
# mxPluginExec/plugins/EDPluginGroupCCP4-v1.0/plugins/EDPluginExecPointlessv1_0.py

import os
import re
from pathlib import Path

from edna2.tasks.AbstractTask import AbstractTask

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging
import traceback
import subprocess

logger = UtilsLogging.getLogger()


class AimlessTask(AbstractTask):
    """
    Execution of CCP4 aimless
    """

    def run(self, inData):
        outData = {}
        input_file = inData['input_file']
        output_file = inData['output_file']
        symoplib = UtilsConfig.get('CCP4', 'symoplib')
        ccp4setup = UtilsConfig.get('CCP4', 'ccp4setup')
        if ccp4setup is None:
            commandLine = ""
        else:
            commandLine = ". " + ccp4setup + '\n'
        if symoplib is None:
            commandLine += 'aimless HKLIN {0} HKLOUT {1}'.format(
                input_file, output_file)
        else:
            commandLine += 'aimless HKLIN {0} HKLOUT {1} SYMINFO {2}'.format(
                input_file, output_file, symoplib)
        logger.info("Command line: {0}".format(commandLine))
        start_image = inData['start_image']
        end_image = inData['end_image']
        projectName = inData.get('dataCollectionID', 'EDNA_proc')
        resolution = inData.get('res', 0.0)
        anom = inData['anom']
        listCommand = [
            'bins 15',
            'run 1 batch {0} to {1}'.format(start_image, end_image),
            'name run 1 project {0} crystal DEFAULT dataset NATIVE'.format(projectName),
            'scales constant',
            'resolution 50 {0}'.format(resolution),
            'cycles 100',
            'anomalous {0}'.format('ON' if anom else 'OFF'),
            'output MERGED UNMERGED',
            'END'
            ]
        self.setLogFileName('aimless.log')
        self.runCommandLine(commandLine, listCommand=listCommand)
        outData['isSuccess'] = os.path.exists(output_file)

        aimlessMergedMtz = self.getWorkingDirectory() / (output_file)
        aimlessUnmergedMtz = self.getWorkingDirectory() / (output_file.replace('.mtz' , "_unmerged.mtz"))
        aimlessLog = self.getWorkingDirectory() / self.getLogFileName()
        #gzip the unmerged aimless.mtz file
        try:
            logger.debug("gzip'ing aimless unmerged file {0}".format(str(aimlessUnmergedMtz)))
            subprocess.call(['gzip', str(aimlessUnmergedMtz)])
        except Exception:
            logger.debug("gzip'ing the file failed: {0}".format(traceback.format_exc()))
        aimlessUnmergedMtzGz = str(aimlessUnmergedMtz) + ".gz"
        outData["aimlessMergedMtz"] = str(aimlessMergedMtz)
        outData["aimlessUnmergedMtz"] = aimlessUnmergedMtzGz
        outData["aimlessLog"] = aimlessLog
        outData['aimlessResults'] = self.extractAimlessResults(aimlessLog)
        return outData
    
    @staticmethod
    def extractAimlessResults(logfile):
        """
        extract the aimless results summary
        """
        aimlessResults = {
            "overall" : {
                "AutoProcScalingStatisticsId" : None,
            },
            "innerShell": {
                "AutoProcScalingStatisticsId" : None,
            },
            "outerShell": {
                "AutoProcScalingStatisticsId" : None,
            },
        }
        extract = []
        try:
            with open(logfile,"r") as fp:
                for line in fp:
                    if line.startswith("<!--SUMMARY_BEGIN--> $TEXT:Result: $$ $$"):
                        while not line.startswith("$$ <!--SUMMARY_END-->"):
                            extract.append(line.strip('\n'))
                            line = next(fp)
                        break
        except:
            logger.error("aimless log file could not be parsed")
            return None
        
        lowResLimit = [x for x in extract if x.startswith("Low resolution limit")][0].split()[-3:]
        aimlessResults["overall"]["resolutionLimitLow"] = float(lowResLimit[0])
        aimlessResults["innerShell"]["resolutionLimitLow"] = float(lowResLimit[1])
        aimlessResults["outerShell"]["resolutionLimitLow"] = float(lowResLimit[2])
        hiResLimit = [x for x in extract if x.startswith("High resolution limit")][0].split()[-3:]
        aimlessResults["overall"]["resolutionLimitHigh"] = float(hiResLimit[0])
        aimlessResults["innerShell"]["resolutionLimitHigh"] = float(hiResLimit[1])
        aimlessResults["outerShell"]["resolutionLimitHigh"] = float(hiResLimit[2])
        rMerge = [x for x in extract if x.startswith("Rmerge  (all I+ and I-)")][0].split()[-3:]
        aimlessResults["overall"]["rmerge"] = float(rMerge[0])
        aimlessResults["innerShell"]["rmerge"] = float(rMerge[1])
        aimlessResults["outerShell"]["rmerge"] = float(rMerge[2])
        rmeasWithinIplusIminus = [x for x in extract if x.startswith("Rmeas (within I+/I-)")][0].split()[-3:]
        aimlessResults["overall"]["rmeasWithinIplusIminus"] = float(rmeasWithinIplusIminus[0])
        aimlessResults["innerShell"]["rmeasWithinIplusIminus"] = float(rmeasWithinIplusIminus[1])
        aimlessResults["outerShell"]["rmeasWithinIplusIminus"] = float(rmeasWithinIplusIminus[2])
        rmeasAllIplusIminus = [x for x in extract if x.startswith("Rmeas (all I+ & I-)")][0].split()[-3:]
        aimlessResults["overall"]["rmeasAllIplusIminus"] = float(rmeasAllIplusIminus[0])
        aimlessResults["innerShell"]["rmeasAllIplusIminus"] = float(rmeasAllIplusIminus[1])
        aimlessResults["outerShell"]["rmeasAllIplusIminus"] = float(rmeasAllIplusIminus[2])
        rpimWithinIplusIminus = [x for x in extract if x.startswith("Rpim (within I+/I-)")][0].split()[-3:]
        aimlessResults["overall"]["rpimWithinIplusIminus"] = float(rpimWithinIplusIminus[0])
        aimlessResults["innerShell"]["rpimWithinIplusIminus"] = float(rpimWithinIplusIminus[1])
        aimlessResults["outerShell"]["rpimWithinIplusIminus"] = float(rpimWithinIplusIminus[2])
        rpimAllIplusIminus = [x for x in extract if x.startswith("Rpim (all I+ & I-)")][0].split()[-3:]
        aimlessResults["overall"]["rpimAllIplusIminus"] = float(rpimAllIplusIminus[0])
        aimlessResults["innerShell"]["rpimAllIplusIminus"] = float(rpimAllIplusIminus[1])
        aimlessResults["outerShell"]["rpimAllIplusIminus"] = float(rpimAllIplusIminus[2])
        nTotalObservations = [x for x in extract if x.startswith("Total number of observations")][0].split()[-3:]
        aimlessResults["overall"]["nTotalObservations"] = int(nTotalObservations[0])
        aimlessResults["innerShell"]["nTotalObservations"] = int(nTotalObservations[1])
        aimlessResults["outerShell"]["nTotalObservations"] = int(nTotalObservations[2])
        nTotalUniqueObservations = [x for x in extract if x.startswith("Total number unique")][0].split()[-3:]
        aimlessResults["overall"]["nTotalUniqueObservations"] = int(nTotalUniqueObservations[0])
        aimlessResults["innerShell"]["nTotalUniqueObservations"] = int(nTotalUniqueObservations[1])
        aimlessResults["outerShell"]["nTotalUniqueObservations"] = int(nTotalUniqueObservations[2])
        meanIoverSigI = [x for x in extract if x.startswith("Mean((I)/sd(I))")][0].split()[-3:]
        aimlessResults["overall"]["meanIoverSigI"] = float(meanIoverSigI[0])
        aimlessResults["innerShell"]["meanIoverSigI"] = float(meanIoverSigI[1])
        aimlessResults["outerShell"]["meanIoverSigI"] = float(meanIoverSigI[2])
        ccHalf = [x for x in extract if x.startswith("Mn(I) half-set correlation CC(1/2)")][0].split()[-3:]
        aimlessResults["overall"]["ccHalf"] = float(ccHalf[0])
        aimlessResults["innerShell"]["ccHalf"] = float(ccHalf[1])
        aimlessResults["outerShell"]["ccHalf"] = float(ccHalf[2])
        completeness = [x for x in extract if x.startswith("Completeness                   ")][0].split()[-3:]
        aimlessResults["overall"]["completeness"] = float(completeness[0])
        aimlessResults["innerShell"]["completeness"] = float(completeness[1])
        aimlessResults["outerShell"]["completeness"] = float(completeness[2])
        multiplicity = [x for x in extract if x.startswith("Multiplicity                   ")][0].split()[-3:]
        aimlessResults["overall"]["multiplicity"] = float(multiplicity[0])
        aimlessResults["innerShell"]["multiplicity"] = float(multiplicity[1])
        aimlessResults["outerShell"]["multiplicity"] = float(multiplicity[2])
        anomalousCompleteness = [x for x in extract if x.startswith("Anomalous completeness")][0].split()[-3:]
        aimlessResults["overall"]["anomalousCompleteness"] = float(anomalousCompleteness[0])
        aimlessResults["innerShell"]["anomalousCompleteness"] = float(anomalousCompleteness[1])
        aimlessResults["outerShell"]["anomalousCompleteness"] = float(anomalousCompleteness[2])
        anomalousMultiplicity = [x for x in extract if x.startswith("Anomalous multiplicity")][0].split()[-3:]
        aimlessResults["overall"]["anomalousMultiplicity"] = float(anomalousMultiplicity[0])
        aimlessResults["innerShell"]["anomalousMultiplicity"] = float(anomalousMultiplicity[1])
        aimlessResults["outerShell"]["anomalousMultiplicity"] = float(anomalousMultiplicity[2])
        ccAno = [x for x in extract if x.startswith("DelAnom correlation between half-sets")][0].split()[-3:]
        aimlessResults["overall"]["ccAno"] = float(ccAno[0])
        aimlessResults["innerShell"]["ccAno"] = float(ccAno[1])
        aimlessResults["outerShell"]["ccAno"] = float(ccAno[2])
        sigAno = [x for x in extract if x.startswith("Mid-Slope of Anom Normal Probability")][0].split()[-3:]
        aimlessResults["overall"]["sigAno"] = float(sigAno[0])
        aimlessResults["innerShell"]["sigAno"] = None
        aimlessResults["outerShell"]["sigAno"] = None

        return aimlessResults

class PointlessTask(AbstractTask):
    """
    Executes the CCP4 program pointless
    """

    def run(self, inData):
        symoplib = UtilsConfig.get('CCP4', 'symoplib')
        ccp4setup = UtilsConfig.get('CCP4', 'ccp4setup')
        logger.debug(f'CCP4 Setup: {ccp4setup}')
        if ccp4setup is None:
            logger.warning('CCP4 setup not found!')
            commandLine = ""
        else:
            commandLine = ". " + ccp4setup + '\n'

        self.input_file = inData['input_file']
        self.output_file = inData['output_file']
        commandLine += 'pointless'
        if UtilsConfig.isEMBL():
            commandLine += ' -c'
        commandLine += " xdsin {0} hklout {1}".format(self.input_file, self.output_file)
        listCommand = ['setting symmetry-based']
        if 'choose_spacegroup' in inData:
            listCommand += 'choose spacegroup {0}'.format(
                inData['choose_spacegroup'])
        self.setLogFileName('pointless.log')
        self.runCommandLine(commandLine, listCommand=listCommand)
        outData = self.parsePointlessOutput(self.getLogPath())
        outData["pointlessUnmergedMtz"] = str(self.getWorkingDirectory() / self.output_file)

        return outData

    @classmethod
    def parsePointlessOutput(cls, logPath):
        sgre = re.compile(""" \* Space group = '(?P<sgstr>.*)' \(number\s+(?P<sgnumber>\d+)\)""")
        outData = {'isSuccess': False}
        if logPath.exists():
            with open(str(logPath)) as f:
                log = f.read()
            m = sgre.search(log)
            if m is not None:
                d = m.groupdict()
                sgnumber = d['sgnumber']
                sgstr = d['sgstr']

                outData['sgnumber'] = int(sgnumber)
                outData['sgstr'] = sgstr
                outData['isSuccess'] = True
                # Search first for unit cell after the Laue group...
                unitCellRe = re.compile("""  Laue group confidence.+\\n\\n\s+Unit cell:(.+)""")
                m2 = unitCellRe.search(log)
                if m2 is None:
                    # Then search it from the end...
                    unitCellRe = re.compile(""" \* Cell Dimensions : \(obsolete \- refer to dataset cell dimensions above\)\\n\\n(.+)""")
                    m2 = unitCellRe.search(log)
                if m2 is not None:
                    listCell = m2.groups()[0].split()
                    cell = {
                        'length_a': float(listCell[0]),
                        'length_b': float(listCell[1]),
                        'length_c': float(listCell[2]),
                        'angle_alpha': float(listCell[3]),
                        'angle_beta': float(listCell[4]),
                        'angle_gamma': float(listCell[5])
                    }
                    outData['cell'] = cell
        return outData

    def gzipUnmergedPointlessFile(self):
        pointless_out = self.getWorkingDirectory() / self.output_file
        try:
            self.DEBUG("gzip'ing pointless multirecord file {0}".format(pointless_out))
            subprocess.call(['gzip', pointless_out])
        except Exception:
            self.DEBUG("gzip'ing the file failed: {0}".format(traceback.format_exc()))

class TruncateTask(AbstractTask):
    """run the CCP4 program Truncate"""

    def getInDataSchema(self):
        return {
             "type": "object",
             "properties": {
                "inputFile": {"type":"string"},
                "outputFile":{"type":"string"},
                "res": {"type":"number"},
                "isAnom": {"type":"boolean"},
                "nres": {"type":["integer","null"]}
             }
        }

    def run(self, inData):
        outData = {}
        ccp4setup = UtilsConfig.get('CCP4', 'ccp4setup')
        logger.debug(f'CCP4 Setup: {ccp4setup}')
        if ccp4setup is None:
            logger.warning('CCP4 setup not found!')
            commandLine = ""
        else:
            commandLine = ". " + ccp4setup + '\n'

        self.inputFile = inData['inputFile']
        self.outputFile = self.getWorkingDirectory() / inData['outputFile']
        commandLine += 'truncate '
        commandLine += 'hklin {0} hklout {1}'.format(self.inputFile, self.outputFile)
        listCommand = ['truncate YES']
        listCommand.append("nres {0}".format(inData.get("nres")) if inData.get("nres") else "")
        listCommand.append('anomalous {0}'.format(inData["isAnom"]))
        listCommand.append('plot OFF')
        listCommand.append('labout F=F_xdsproc SIGF=SIGF_xdsproc')
        listCommand.append('falloff YES')
        listCommand.append('resolution 50 {0}'.format(inData["res"]))
        listCommand.append('PNAME EDNA2proc')
        listCommand.append('DNAME EDNA2proc')
        listCommand.append('end')

        self.setLogFileName('truncate.log')
        logger.debug("Running ccp4/truncate...")
        self.runCommandLine(commandLine, listCommand=listCommand)

        outData["truncateOutputMtz"] = self.outputFile
        outData["truncateLogPath"] = self.getWorkingDirectory() / self.getLogFileName()
        self.isSuccess = Path(self.outputFile).exists()

        return outData
    
class UniqueifyTask(AbstractTask):
    def run(self, inData):
        outData = {}
        ccp4setup = UtilsConfig.get('CCP4', 'ccp4setup')
        logger.debug(f'CCP4 Setup: {ccp4setup}')
        if ccp4setup is None:
            logger.warning('CCP4 setup not found!')
            commandLine = ""
        else:
            commandLine = ". " + ccp4setup + '\n'

        self.inputFile = inData['inputFile']
        self.outputFile = self.getWorkingDirectory() / inData['outputFile']
        self.setLogFileName('uniqeuify.log')

        commandLine += 'uniqueify '
        commandLine += '{0} {1}'.format(self.inputFile, self.outputFile)

        logger.debug("Running ccp4/uniqueify...")

        self.runCommandLine(commandLine)
        outData["uniqueifyOutputMtz"] = self.outputFile
        self.isSuccess = Path(self.outputFile).exists()

        return outData

class DimpleTask(AbstractTask):
    """
    This task runs dimple to replace B-factors and (optinally) break the model into domains
    """

    def run(self, inData):
        output_Dir = self._workingDirectory 
        outData = {}
        if os.environ.get('CCP4', None) is None:
            commandLine = 'source /mxn/groups/sw/mxsw/env_setup/ccp4_env.sh \n'
        else:
            commandLine = ''
            logger.info(f"CCP4 version is {os.environ.get('CCP4_VERSION', None)}")

        commandLine += 'dimple '
        commandLine += f'{inData["PDB_file"]} '
        commandLine += f'{inData["MTZ_file"]} '
        commandLine += f'{output_Dir}'

        logPath = self.getWorkingDirectory() / 'dimple.log'
        self.runCommandLine(commandLine, logPath=logPath)
        with open(str(logPath)) as f:
            logText = f.read()
        
        logger.info("Command line: {0}".format(commandLine))

        # outData = self.parseProcessPredictedModel(logPath)
        
        if Path(f"{output_Dir}/final.pdb").exists() and Path(f"{output_Dir}/final.mtz").exists():
            outData["isSuccess"] = True

        return outData
    
    def parseDimpleLog(self, logPath):
        if logPath.exists():
            with open(str(logPath)) as f:
                log = f.read()

        


