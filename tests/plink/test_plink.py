from copy import deepcopy
import os

from numpy import array
from numpy import nan
import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal
import pytest

import cdpybio as cpb

def add_root(fn):
    return os.path.join(cpb._root, 'tests', 'star', fn)

LINEAR = add_root('test.glm.linear')
LINEAR_REWRITE = add_root('test.glm.linear.rewrite')
LOGISTIC = add_root('test.glm.logistic.hybrid')
LOGISTIC_REWRITE = add_root('test.glm.logistic.hybrid.rewrite')

class TestParsing:
    def test_linear(self):
        vals = array([['1', 86028, '0', '0', 'ADD', 0, nan, nan, nan, nan],
            ['1', 86028, '0', '0', 'age', 0, nan, nan, nan, nan],
            ['1', 86028, '0', '0', 'sex', 0, nan, nan, nan, nan],
            ['1', 86028, '0', '0', 'C1', 0, nan, nan, nan, nan],
            ['1', 86028, '0', '0', 'C2', 0, nan, nan, nan, nan],
            ['1', 86028, '0', '0', 'C3', 0, nan, nan, nan, nan],
            ['1', 86028, '0', '0', 'C4', 0, nan, nan, nan, nan],
            ['1', 723307, '1', '2', 'ADD', 423, -0.067378, 0.191871, -0.351163,
             0.725645],
            ['1', 723307, '1', '2', 'age', 423, -0.394825, 0.046971399999999996,
             -8.40565, 6.875239999999999e-16],
            ['1', 723307, '1', '2', 'sex', 423, 0.38931, 0.045982499999999996,
             8.46648, 4.4102499999999995e-16],
            ['1', 723307, '1', '2', 'C1', 423, -0.043173800000000005,
             0.047640800000000004, -0.906236, 0.36533699999999997],
            ['1', 723307, '1', '2', 'C2', 423, -0.0949844, 0.0449457, -2.11331,
             0.0351702],
            ['1', 723307, '1', '2', 'C3', 423, 0.0532012, 0.04523219999999999,
             1.17618, 0.240198],
            ['1', 723307, '1', '2', 'C4', 423, -0.08434960000000001, 0.0461369,
             -1.82825, 0.0682302]], dtype=object)
        ind = [u'Affx-14150122', u'Affx-14150122', u'Affx-14150122',
               u'Affx-14150122', u'Affx-14150122', u'Affx-14150122',
               u'Affx-14150122', u'Affx-13546538', u'Affx-13546538',
               u'Affx-13546538', u'Affx-13546538', u'Affx-13546538',
               u'Affx-13546538', u'Affx-13546538']
        cols = [u'CHROM', u'POS', u'REF', u'ALT1', u'TEST', u'OBS_CT', u'BETA',
                u'SE', u'T_STAT', u'P']
        df = pd.DataFrame(vals, index=ind, columns=cols)
        df2 = cpb.plink.read_linear2(LINEAR)
        assert_frame_equal(df, df2)
