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

__authors__ = ['O. Svensson']
__license__ = 'MIT'
__date__ = '21/04/2019'

from edna2.utils import UtilsLogging

from edna2.tasks.AbstractTask import AbstractTask

from tests.test_tasks.ParallelTasks.ExampleTask import ExampleTask

logger = UtilsLogging.getLogger()


class ControlTestTask(AbstractTask):
    """
    Control task for testing parallel execution
    """

    def run(self, inData):
        numberOfTasks = inData['numberOfTasks']
        startNumber = inData['startNumber']
        listTasks = []
        for taskNumber in range(startNumber, startNumber+numberOfTasks):
            testTaskInData = {'taskNumber': taskNumber}
            testTask = ExampleTask(inData=testTaskInData)
            logger.info(testTask)
            listTasks.append(testTask)
            testTask.start()
        for testTask in listTasks:
            testTask.join()
        outData = {'status': 'finished'}
        return outData
