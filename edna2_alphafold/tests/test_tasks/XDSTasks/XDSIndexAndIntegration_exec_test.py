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
__date__ = "20/04/2020"

import os
import unittest

from edna2.utils import UtilsTest
from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

from edna2.tasks.XDSTasks import XDSIndexAndIntegration

logger = UtilsLogging.getLogger()


class XDSItegrationExecTest(unittest.TestCase):

    def setUp(self):
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)

    # @unittest.skipIf(UtilsConfig.getSite() == 'Default',
    #                  'Cannot run test with default config')
    @unittest.skipIf(True, "Temporarily disabled due to missing test data")
    def test_execute_XDSIntegration(self):
        referenceDataPath = self.dataPath / 'id30a1_1_fast_char.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        xdsIntegration = XDSIndexAndIntegration(inData=inData)
        xdsIntegration.execute()
        self.assertTrue(xdsIntegration.isSuccess())
