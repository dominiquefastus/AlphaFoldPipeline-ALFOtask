#
# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the 'Software'), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__authors__ = ["D. Fastus"]
__license__ = "MIT"
__date__ = "16/12/2023"

import gemmi

from edna2.utils import UtilsConfig
from edna2.utils import UtilsLogging

def construct_ensembles(json_inputs):
    logger = UtilsLogging.getLogger()
    logger.info('Starting to construct ensembles')

    title = json_inputs['title']
    mode = json_inputs['mode']
    root = json_inputs['root']
    mtz_file = json_inputs['mtz_file']
    ensembles = json_inputs['ensembles']

    logger.info(f'Reading .mtz file: {mtz_file}')
    mtz = gemmi.read_mtz_file(mtz_file)
    f_col = next((col.label for col in mtz.columns if col.type == 'F'), None)
    sigf_col = next((col.label for col in mtz.columns if col.type == 'Q'), None)
    i_col = next((col.label for col in mtz.columns if col.type == 'J'), None)
    sigi_col = next((col.label for col in mtz.columns if col.type == 'L'), None)

    logger.info('Determining LABIn line based on available columns')
    if i_col and sigi_col:
        labin_line = f'LABIn I={i_col} SIGI={sigi_col}'
    elif f_col and sigf_col:
        labin_line = f'LABIn F={f_col} SIGF={sigf_col}'
    else:
        logger.error('No suitable columns found in .mtz file')
        raise ValueError('No suitable columns found in .mtz file')

    logger.info('Constructing Phaser script')
    phaser_script = f'''phaser << eof
TITLE {title}
MODE {mode}
HKLIN {mtz_file}
{labin_line}
'''

    ensemble_lines = ''
    composition_lines = ''
    search_lines = ''

    for ensemble_info in ensembles:
        ensemble = ensemble_info['ensemble']
        pdb = ensemble_info['pdb']
        identity = ensemble_info['identity']
        sequence = ensemble_info['sequence']
        num = ensemble_info['num']

        ensemble_lines += f'ENSEMBLE {ensemble} PDBFILE {pdb} IDENTITY {identity}\n'
        composition_lines += f'COMPOSITION PROTEIN SEQUENCE {sequence} NUM {num}\n'
        search_lines += f'SEARCH ENSEMBLE {ensemble} NUM {num}\n'

    phaser_script += ensemble_lines + composition_lines + search_lines
    phaser_script += f'ROOT {root} # not the default\neof'

    logger.info('Finished constructing Phaser script')
    return phaser_script