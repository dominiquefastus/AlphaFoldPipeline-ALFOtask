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

import pprint
import unittest
import tempfile

from edna2.utils import UtilsTest
from edna2.utils import UtilsLogging

from edna2.tasks.XDSTasks import XDSTask
from edna2.tasks.XDSTasks import XDSIndexing

logger = UtilsLogging.getLogger()


class XDSTasksUnitTest(unittest.TestCase):

    def setUp(self):
        self.dataPath = UtilsTest.prepareTestDataPath(__file__)

    def tes_writeSPOT_XDS(self):
        spotFile = self.dataPath / '00001.spot'
        spotXdsReferenceFile = self.dataPath / 'SPOT.XDS'
        with open(str(spotXdsReferenceFile)) as f:
            spotXdsReference = f.read()
        spotXds = XDSTask.createSPOT_XDS([spotFile], oscRange=1)
        with open('/tmp/SPOT.XDS', 'w') as f:
            f.write(spotXds)
        with open('/tmp/SPOT.XDS_REF', 'w') as f:
            f.write(spotXdsReference)
        self.assertEqual(spotXdsReference.split('\n')[0], spotXds.split('\n')[0])
        self.assertEqual(spotXdsReference, spotXds)

    def test_readIdxrefLp(self):
        idxRefLpPath = self.dataPath / 'IDXREF.LP_TRYP'
        resultXDSIndexing = XDSIndexing.readIdxrefLp(idxRefLpPath)

    def test_parseXparm(self):
        xparmPath = self.dataPath / 'XPARM.XDS'
        xparmDict = XDSIndexing.parseXparm(xparmPath)
        self.assertIsNotNone(xparmDict)

    def test_getXDSDetector(self):
        referenceDataPath = self.dataPath / 'inDataXDSIndexing.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        firstSubWedge = inData["subWedge"][0]
        dictDetector = firstSubWedge["experimentalCondition"]["detector"]
        dictXDSDetector = XDSTask.getXDSDetector(dictDetector)
        # pprint.pprint(dictXDSDetector)
        self.assertTrue(dictXDSDetector["name"] == "PILATUS")

    def test_generateXDS_INP(self):
        referenceDataPath = self.dataPath / 'inDataXDSIndexing.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        listXDS_INP = XDSTask.generateXDS_INP(inData)
        pprint.pprint(listXDS_INP)

    def test_generateImageLinks_1(self):
        referenceDataPath = self.dataPath / 'inDataXDSIntegration.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        XDSTask.generateImageLinks(inData)

    def test_generateImageLinks_2(self):
        referenceDataPath = self.dataPath / 'inDataXDSIntegration_one_subWedge.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        XDSTask.generateImageLinks(inData)

    def test_generateImageLinks_3(self):
        referenceDataPath = self.dataPath / 'inDataXDSGenerateBackground_eiger16m.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        imageLinks = XDSTask.generateImageLinks(inData)
        pprint.pprint(imageLinks)

    def test_generateImageLinks_3(self):
        referenceDataPath = self.dataPath / 'id30a1_1_fast_char.json'
        inData = UtilsTest.loadAndSubstitueTestData(referenceDataPath)
        imageLinks = XDSTask.generateImageLinks(inData)
        pprint.pprint(imageLinks)
