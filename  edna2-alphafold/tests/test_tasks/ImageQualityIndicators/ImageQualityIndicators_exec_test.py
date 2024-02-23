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

import os
import unittest

from edna2.utils import UtilsTest
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

from edna2.tasks.ImageQualityIndicators import ImageQualityIndicators

logger = UtilsLogging.getLogger()


class ImageQualityIndicatorsExecTest(unittest.TestCase):
    def setUp(self):
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)
        # self.dataPath = pathlib.Path(os.getcwd()) / 'data'

    @unittest.skipIf(
        UtilsConfig.getSite() == "Default",
        "Cannot run ImageQualityIndicatorsExecTest " + "test with default config",
    )
    @unittest.skipIf(
        not os.path.exists(
            "/data/scisoft/pxsoft/data/WORKFLOW_TEST_DATA/id30a2/inhouse/opid30a2"
            + "/20191129/RAW_DATA/t1/MeshScan_05/mesh-t1_1_0001.cbf"
        ),
        "Cannot find CBF file mesh-t1_1_0001.cbf",
    )
    def test_execute(self):
        referenceDataPath = self.dataPath / "inDataImageQualityIndicatorsTask.json"
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        task = ImageQualityIndicators(inData=inData)
        task.execute()
        self.assertFalse(task.isFailure())
        outData = task.outData
        self.assertTrue("imageQualityIndicators" in outData)
        # self.assertTrue('resolution_limit' in outData['crystfel_results'][0])
        self.assertEqual(72, len(outData["imageQualityIndicators"]))


if __name__ == "__main__":
    unittest.main()
