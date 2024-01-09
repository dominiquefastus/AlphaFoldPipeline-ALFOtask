import unittest
import json
from UtilsPhaserScript import construct_ensembles

class UtilsPhaserScriptTest(unittest.TestCase):
    def test_construct_ensembles(self):
        with open('inDataPhenixPhaser.json', 'r') as f:
            json_inputs = json.load(f)

        expected_output = '''phaser << eof
TITLe Test Title
MODE Test Mode
HKLIn Test File
LABIn F=None SIGF=None
ENSEmble Test Ensemble PDB Test PDB IDENtity Test Identity
COMPosition PROTein SEQuence Test Sequence NUM Test Num #Test Ensemble
SEARch ENSEmble Test Ensemble NUM Test Num
ROOT Test Root # not the default\neof'''
        result = construct_ensembles(json_inputs)
        self.assertEqual(result, expected_output)

if __name__ == '__main__':
    unittest.main()